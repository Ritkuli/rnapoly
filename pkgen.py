#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" Generates a list of kissing loop pairs. Also contains methods to read
a generated file and select orthogonal kissing loops.

Author: Antti Elonen
"""
from tempfile import NamedTemporaryFile
import itertools, threading, traceback, shutil, re
import  argparse, os, math, sys, random, subprocess


PFUNC = shutil.which("pfunc")
if not PFUNC: PFUNC = os.path.join(script_path, "lib/nupack/pfunc")
CMD = (PFUNC,)
THREADS = 12
TEMP = ".temp.txt"
CHUNK_SIZE = 1024
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "W": "AU",
    "S": "CG",
    "M": "AC",
    "K": "GU",
    "R": "AG",
    "Y": "CU",
    "B": "CGU",
    "D": "AGU",
    "H": "ACU",
    "V": "ACG",
    "N": "ACGU",
}


def main(argv):
    args = interpret_arguments(argv)
    print(args.task, args.output)
    if args.task.lower() == "generate":
        strings = generate_strings(args.length)
        pairs = get_pairs(strings, args.count)
        output = save(pairs, args.output)
    elif args.task.lower() == "select":
        choose_pairs(args.pairs_file, method = args.method, count = args.count)
    else:
        print("Invalid task.")
    return


def interpret_arguments(argv):
    """Interprets the commandline arguments

    Returns:
        namespace object
    """

    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(dest="task", required=True)

    parser_gen = sub_parsers.add_parser("generate")
    parser_gen.add_argument("output", type=str,
                    help="The output path.")
    parser_gen.add_argument("-l", "--length", default = 6, type=int,
                    help="Pseudoknot length")
    parser_gen.add_argument("-c", "--count", default = 8, type=int,
                    help="Number of pseudoknots to generate / select.")

    parser_select = sub_parsers.add_parser("select")
    parser_select.add_argument("pairs_file", type=str,
                    help="The pairs file from which to select pseudoknots.")
    parser_select.add_argument("method", type=str,
                    help="The method used to select the pseudoknots.", choices=("anneal", "list"))
    parser_select.add_argument("output", type=str,
                    help="The output path.")
    parser_select.add_argument("-l", "--length", default = 6, type=int,
                    help="Pseudoknot length")
    parser_select.add_argument("-c", "--count", default = 8, type=int,
                    help="Number of pseudoknots to generate / select.")

    args = parser.parse_args(argv)
    return args

def in_conflict(strand, prevent, with_180 = True):
    """ Checks whether the given strand is in conflict with the prevent sequences.

    Args:
        strand -- str: IUPAC primary structure
        prevent -- list<str>: prevent these sequences
    KWArgs:
        with_180 -- bool, with UGAA + ACG bases around the strand
    Returns:
        bool
    """
    for p in prevent:
        m = re.compile(p)
        if with_180: s = "UGAA" + strand + "ACG"
        else: s = strand
        if m.search(s):
            #print(s, m.search(s), m.search(s).re)
            return True
    return False

def get_prevent_regex(prevent_str):
    """ Returns a regular expression string to match the prevent sequences.

    Args:
        prevent_str -- str: regular expression
    """
    prevent = []
    if prevent_str == None: return prevent
    for p in prevent_str.split(","):
        p = p.strip()
        current = []
        for c in p:
            t = "(" + "|".join(list(IUPAC[c])) + ")"
            current.append(t)
        prevent.append("".join(current))
    #print(prevent)
    return prevent


def choose_pairs(pairs_f, threshold = 3.75, count = 8, GC_content = 1.0, prevent = None, enforced = [], method = "list"):
    """ Returns kissing loop pairs selected from a kissing loops file according to the
    input parameters.

    Args:
        pairs_f -- str: path to the kissing loops pairs file
    Kwargs:
        threshold -- float: energy threshold between two pairs to consider them orthogonal
        count -- int: number of pairs to select
        GC_content -- float: maximum GC_content in the pairs
        prevent -- list<str>: prevent from choosing anything matching to these
        enforced -- list<tuple<str, str>>: always select these pairs
        method -- str: the method used to select the kissing loops
    Returns:
        list<str> -- kissing loops
    """
    if method == "list":
        return list_choose(pairs_f, threshold, count, GC_content, prevent, enforced)
    elif method == "anneal":
        return anneal_choose(pairs_f, threshold, count, GC_content, prevent, enforced)
    else:
        from inspect import currentframe, getframeinfo
        frameinfo = getframeinfo(currentframe())
        raise(Exception("Invalid method in {} {}".format(frameinfo.filename, frameinfo.lineno)))

def anneal_choose(pairs_f, threshold, count, GC_content, prevent, enforced, dT = 0.99, step = 20, max_energy = -7):
    """ Chooses n kissing loops based on the input parameters using simulated annealing.

    Args:
        pairs_f -- str: path to the kissing loops pairs file
    Kwargs:
        threshold -- float: energy threshold between two pairs to consider them orthogonal
        count -- int: number of pairs to select
        GC_content -- float: maximum GC_content in the pairs
        prevent -- list<str>: prevent from choosing anything matching to these
        enforced -- list<tuple<str, str>>: always select these pairs
        method -- str: the method used to select the kissing loops
    Returns:
        list<str> -- kissing loops
    """
    prevent_regex = get_prevent_regex(prevent)
    pfg = read_pair_file(pairs_f, GC_content)
    d = {}
    print("Reading pairs file...")
    max_i = 0
    for i, (pair, val) in enumerate(pfg):
        if max_i == 0 and val > max_energy: max_i = i
        elif val > max_energy * 0.5: break
        d[pair] = val
    keys = list(d.keys())
    pairs = keys[:count]
    print("Annealing...")
    T = 1000
    trajectory = []
    best = [None, float("inf")]
    i = 0
    while T > 0.1:
        if not i % step:
            T = dT * T
        i += 1
        # get candidate
        drop_i = random.randint(0, len(pairs) - 1)
        c_pair = get_candidate(keys, prevent_regex, max_i)
        c_pairs = pairs[:drop_i] + pairs[drop_i + 1:]
        can_cost = get_cost(c_pair, c_pairs, d)
        cur_cost = get_cost(pairs[drop_i], c_pairs, d)

        # test and implement candidate
        p = math.exp(-(can_cost - cur_cost) / T)
        if(random.random() < p):
            pairs[drop_i] = c_pair

        tot = sum(get_all_costs(pairs, d))
        trajectory.append(tot)
        #print("Total cost: {}. \nPairs: {}\n-----------------".format(tot, pairs))
        if tot < best[1]: best = ([p for p in pairs], tot)

    #print(trajectory)
    assert len(best[0]) == count
    costs = get_all_costs(best[0], d)
    print("Score:", sum(costs))
    print("Costs:", costs)
    for pair in best[0]:
        print("Chosen: {} {}".format(*pair))
    return best[0]


def get_candidate(keys, prevent_regex, max_i):
    """

    Args:

    Returns:
    """
    while True:
        c_pair = keys[random.randint(0, max_i)]
        if not is_complement(*c_pair) \
                or in_conflict(c_pair[0], prevent_regex) \
                or in_conflict(c_pair[1][::-1], prevent_regex):
            continue
        else:
            return c_pair


def get_all_costs(pairs, d):
    """

    Args:

    Returns:
    """
    costs = []
    for i in range(len(pairs)):
        c_pairs = pairs[:i] + pairs[i + 1:]
        costs.append(get_cost(pairs[i], c_pairs, d))
    return costs


def get_cost(c_pair, pairs, d):
    """

    Args:

    Returns:
    """
    val = d[c_pair]
    mismatch_cost = 0
    for pair in pairs:
        mismatches = itertools.product(c_pair, pair)
        cost = min((d[m] if m in d else float("inf") for m in mismatches))
        if cost < mismatch_cost: mismatch_cost = cost
    c_cost = val - mismatch_cost
    return c_cost


def list_choose(pairs_f, threshold, count, GC_content, prevent, enforced):
    """ Chooses n kissing loops based on the input parameters using list traversal.

    Args:
        pairs_f -- str: path to the kissing loops pairs file
    Kwargs:
        threshold -- float: energy threshold between two pairs to consider them orthogonal
        count -- int: number of pairs to select
        GC_content -- float: maximum GC_content in the pairs
        prevent -- list<str>: prevent from choosing anything matching to these
        enforced -- list<tuple<str, str>>: always select these pairs
        method -- str: the method used to select the kissing loops
    Returns:
        list<str> -- kissing loops
    """
    prevent_regex = get_prevent_regex(prevent)
    chosen = list(enforced[:min(count, len(enforced))])
    used = set()
    pfg = read_pair_file(pairs_f, GC_content)
    start = next(pfg)
    pairs = [start]
    d = {start[0][0]: [(start[0][1], start[1])], start[0][1]: [(start[0][0], start[1])]}
    enforced_found = 0
    vals = {}
    for pair, val in pairs:
        vals[pair] = val
        if len(chosen) >= count or val > 0: break
        p = set(d[pair[0]]).union(d[pair[1]])
        conflict = False
        if [pair[0], pair[1]] in enforced or [pair[1], pair[0]] in enforced:
            enforced_found += 1# Makes sure enforced pairs are already seen
            print("Enforced: {}".format(pair))
        if enforced_found < len(enforced): conflict = True
        #print(pair[0], in_conflict(pair[0], prevent_regex))
        #print(pair[1], in_conflict(pair[1][::-1], prevent_regex))
        if not is_complement(*pair) or in_conflict(pair[0], prevent_regex) or in_conflict(pair[1][::-1], prevent_regex):
            conflict = True
        while pairs[-1][1] - threshold < val:
            pair_t, val_t = next(pfg)
            pairs.append((pair_t, val_t))
            d.setdefault(pair_t[0], []).append((pair_t[1], val_t))
            d.setdefault(pair_t[1], []).append((pair_t[0], val_t))
        for x in p:
            if conflict: break
            if x[0] in used and val > x[1] - threshold:
                conflict = True
        if not conflict:
            chosen.append(pair)
            used.add(pair[0])
            used.add(pair[1])
            print("Found:", pair, " ::: ", val)
    #print(chosen)
    assert len(chosen) == count
    # print("Score:", get_all_costs(chosen, vals)) # not a good representation without all the pairs
    for pair in chosen:
        print("Chosen: {} {}".format(*pair))
    return chosen

def read_pair_file(pairs_f, GC_content = 1.0):
    """ Reads a pairs file and returns a kissing loop generator

    args:
        pairs_f -- str: path to the pairs
    KWargs:
        GC_content -- maximum proportion of U's and G's in the kissing loop to yield
    returns:
        str -- space separated kissing loop pair
        float -- the energy of the pair
    """
    with open(pairs_f) as f:
        for l in f:
            t = l.strip().split(" ")
            pair = tuple(t[:-1])
            s = "".join(pair)
            frac = (s.count("G") + s.count("C")) / len(s)
            if frac >= GC_content: continue
            val = float(t[-1])
            yield (pair, val)

def get_pairs(strings, count = float("inf")):
    """ Calculates the energies for all possible pairs with the given loops

    args:
        strings -- list<str>: a list of loops
    KWargs:
        count -- int: the number of kissing loop pairs to generate
    returns:
        list<tuple<str, float>> -- a sorted list of the kissing loop pairs and energies
    """
    pairs = []
    threads = THREADS * [None]
    stack = []
    t_file = open(TEMP, "w+")
    processed = 0
    for pair in product(strings, strings):
        if processed > count: break
        processed += 1
        stack.append(pair)
        if len(stack) < THREADS: continue
        print("{} / {} : {} %".format(processed, count, 100 * processed / count))
        for i in range(THREADS):
            p = stack.pop()
            t = threading.Thread(target=get_fitness_and_save, args=(pairs, p))
            threads[i] = t
            t.start()
        for t in threads:
            t.join()
        if not processed % CHUNK_SIZE:
            flush(t_file, pairs)
            pairs.clear()
    flush(t_file, pairs)
    return get_sorted(t_file)

def flush(file, pairs):
    """ Flushes the given pairs to an output file

    args:
        file -- str: path to the output file
        pairs: list<tuple<str, str, float>>: kissing loops and their energies
    """
    for p in pairs:
        file.write("{} {} {}\n".format(*p[0], p[1]))


def get_sorted(file):
    """ Sorts the kissing loop pairs in a file

    args:
        file -- str: path to the file with the kissing loops
    returns:
        list<tuple<str, str, float>> -- sorted list of the kissing loop pairs based on their energies
    """
    file.seek(0)
    pairs = []
    for l in file:
        print(l)
        t = l.strip().split(" ")
        pairs.append(((t[:2]), float(t[2])))
    return sorted(pairs, key = lambda x: x[1])

def get_fitness_and_save(container, pair):
    """ Calculates the fitness, or energy, of the given pair and saves it in the container

    args:
        container -- list-like
        pair -- tuple<str, str>
    """
    v = get_fitness(*pair)
    container.append((pair, v))

def save(pairs, path):
    """ Saves thbest[0]e pairs in the output file

    args:
        pairs -- list<tuple<str, str, float>>
        path -- str: path to the output file
    returns:
        str -- the path of the output file
    """
    path = path if "".endswith(".txt") else path + ".txt"
    with open(path, mode='w') as f:
        for p in pairs:
            f.write("{} {} {}\n".format(*p[0], p[1]))
    print("Successfully saved to {}".format(path))
    return path

def generate_strings(length = 6):
    """ Generates all possible strings of length N using the alphabet AUGC

    args:
        length -- int
    returns:
        list<str>
    """
    combinations = [p for p in itertools.combinations_with_replacement("AUGC", length)]
    strings = set()
    for c in combinations:
        for p in itertools.permutations(c, length):
            strings.add("".join(p))
    return strings

def product(l1, l2):
    """ Returns all possible combinations of entries in two arrays

    Args:
        l1 -- list<str>
        l2 -- list<str>
    Returns:
        iterable
    """
    return itertools.product(l1, l2)


def get_fitness(s1, s2, add_padding = True):
    """ Calcualtes the energy of the given kissing loop pair

    args:
        s1 -- str: first loop
        s2 -- str: second loop
    KWargs:
        add_padding -- bool: add UGAA + AGC padding around the loops
    returns:
        float
    """
    #TODO: add_padding
    t = generate_in(s1, s2)
    if not t: return float("inf")
    name = t.name[:-3]
    cmd = (*CMD, name)
    try:
        a = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output = a.communicate()[0]
        #print(output)
        result = float(re.sub(".*nergy.*:", "", str(output)).split("\\n")[1])
        #print(result)
        return result
    except Exception as e:
        #print(e)
        return float("inf")

def is_base_complement(c1, c2):
    """Returns true if c1 is a watson-crick complement of c2

    args:
        c1 -- char
        c2 -- char
    returns:
        bool -- True if is complement, false otherwise
    """
    if c1 == "A":
        return c2 == "U"
    elif c1 == "U":
        return c2 == "G" or c2 == "A"
    elif c1 == "G":
        return c2 == "C" or c2 == "U"
    elif c1 == "C":
        return c2 == "G"
    return False

def is_complement(s1, s2):
    """Returns true if s1 is a watson-crick complement of s2

    args:
        c1 -- str
        c2 -- str
    returns:
        bool -- True if is complement, false otherwise
    """
    for c1, c2 in zip(s1, s2):
        if not is_base_complement(c1, c2): return False
    return True


def get_structure(s1, s2):
    """Generates a dotbracket structure of two strands

    args:
        strand1 -- str
        strand2 -- str
    returns:
        str -- The dotbracket structure
    """
    result = []
    for i, c in enumerate(s1):
        result.append("(") if is_base_complement(c, s2[i]) else result.append(".")
    complement = [c.replace("(", ")") for c in result[::-1]]
    return "+".join(("".join(result), "".join(complement)))

def generate_in(s1, s2):
    """Generates a temporary in file for NUPACK

    args:
        strand1 -- str
        strand2 -- str
    returns:
        str -- The path to the temporary np-file. None, if flawed input
    """
    #structure = get_structure(s1, s2)
    #if structure.count("(") == 0: return None
    #print("{}, {}: {}".format(s1, s2, structure))
    t = NamedTemporaryFile(suffix = ".in", mode = "w")
    #contents = "2\n{}\n{}\n1 2\n{}".format(s1, s2, structure)
    contents = "{}+{}".format(s1, s2[::-1])
    t.write(contents)
    t.seek(0)
    return t

if __name__ == "__main__":
    main(sys.argv[1:])
