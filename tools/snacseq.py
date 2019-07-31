#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" Generates the primary structure of a SNAC file.

Author: Antti Elonen
"""

from lib.utils import read_snac, write_snac
from tempfile import NamedTemporaryFile
from shutil import copyfile, SameFileError, which
import argparse, sys, subprocess, itertools, re, pkgen, os

script_path = os.path.dirname(sys.argv[0])

OUTPUT_PAT = "output/"
PRIMARY_GENERATOR = which("multitubedesign")
if not PRIMARY_GENERATOR: PRIMARY_GENERATOR = script_path + "/lib/nupack/multitubedesign"
PAIR_LIST = script_path + "/resources/pairs/default.txt"
NP_TEMPLATE = script_path + "/resources/templates/default.np"

def main(argv):
    args = interpret_arguments(argv)
    #print(args.snac)

    d, snac_path = read_snac(args.snac)
    positions = d.get("positions", "")
    orig_primary_structure = d.get("primary_structure", "")
    secondary_structure = d.get("secondary_structure", "")
    pseudoknot_numbering = d.get("pseudoknot_numbering", "")
    length = 6 #TODO: don't treat all pseudoknots as equal lengths
    #print(pseudoknot_numbering)
    count_pseudoknots = secondary_structure.count(length * "[")
    if args.ignore_conflicts:
        prevent = None
    else:
        with open(NP_TEMPLATE) as f:
            for l in f:
                if l.startswith("prevent"):
                    prevent = l.split("=")[-1].strip()
                    break
    #assign pseudoknots and reverse the second pseudoknots
    #TODO: make into a function and make pkgen return the second pseudoknot 5'-3'
    enforced = read_enforced_pairs(args.enforced_pairs)
    pairs = pkgen.choose_pairs(args.pair_list, count = count_pseudoknots, threshold = args.threshold, GC_content = args.gc_content, prevent=prevent, enforced=enforced)
    pseudoknots = []
    visited = {}
    for n in pseudoknot_numbering.split(" "):
        if len(n) < 1: continue
        if n not in visited:
            pair = pairs.pop()
            visited[n] = pair
            pseudoknots.append(pair[0])
        else:
            pseudoknots.append(visited[n][1][::-1])
    #print(pseudoknots)
    np, np_input = generate_np(secondary_structure, pseudoknots, ignore_conflicts=args.ignore_conflicts, restriction=orig_primary_structure, ignore_primary=args.ignore_primary)
    d["nupack_input"] = " | ".join(split_to_domains(np_input, secondary_structure))
    if not args.only_pseudoknots:
        primary_structure = generate_primary(np, args.primary_generator)
        d["primary_structure"] = primary_structure
    if args.output:
        output = args.output
        if not output.endswith(".snac"): output += ".snac"
    elif args.replace:
        output = snac_path
    else:
        output = args.snac.replace(".snac", "") + "_seq.snac"
    try:
        copyfile(snac_path, output)
    except SameFileError: pass
    s = write_snac(output, d)
    print("Saved output to {}".format(s))
    return


def interpret_arguments(argv):
    """Interprets the console arguments
    """
    parser = argparse.ArgumentParser()
    opt = parser.add_argument_group()

    opt.add_argument("snac", type=str,
                     help="The input file in snac format.")
    opt.add_argument("-p", "--pair_list", default=PAIR_LIST, type=str,
                     help="")
    opt.add_argument("-e", "--enforced_pairs", default=None, type=str,
                     help="")
    opt.add_argument("-o", "--output", default=None, type=str,
                     help="The output file name. If not specified, output is input file + _seq.")
    opt.add_argument("-r", "--replace", action="store_true",
                     help="Replace the input file with the output.")
    opt.add_argument("-c", "--ignore_conflicts", action="store_true",
                     help="Ignore conflicts with prevented sequences.")
    opt.add_argument("-i", "--ignore_primary", action="store_true",
                     help="Ignore the existing primary structure. Otherwise use it as a restriction in generation.")
    opt.add_argument("-gc", "--gc_content", default=1.0, type=float,
                     help="The maximum gc content of a pseudoknot in percents.")
    opt.add_argument("-t", "--threshold", default=3.0, type=float,
                     help="The minimum threshold between mismatching pseudoknots. Default = 3.")
    opt.add_argument("-x", "--only_pseudoknots", action="store_true",
                     help="Don't run primary structure design. Only find pseudoknots.")
    opt.add_argument("-P", "--primary_generator", default=PRIMARY_GENERATOR, type=str,
                     help="The executable for the primary structure generator")

    args = parser.parse_args(argv)
    return args

def read_enforced_pairs(path):
    """ Reads a file containing kissing loop pairs into an array.

    Args:
        path -- str
    Returns:
        list<i, j>
    """
    pairs = []
    if path != None:
        with open(path) as f:
            for line in f:
                t = line.split(" ")
                pairs.append(([t[0], t[1]]))
    return pairs


def generate_primary(np, primary_generator):
    """Generates a primary structure sans pseudoknots with the given generator

    Args:
        np -- str: nupack input
        primary_generator -- str: The command to run the primary-
            generating software in commandline.
    Returns:
        str -- The primary structure if successful, None If unsuccessful
    """
    try:
        print("Generating primary structure. This will take a while...")
        t = NamedTemporaryFile(suffix = ".np", mode = "w")
        t.write(np)
        t.seek(0)
        name = t.name[:-3]
        subprocess.call((primary_generator, name))
        output = name + "_0.npo"
        primary = read_primary(output)
        print("Successfully generated primary structure: %s" % primary)
        return primary
    except Exception as e:
        print(e)
        print("Could not run NUPACK tubedesign. Make sure you have NUPACK installed or use the online version")
        sys.exit(0)


def create_domain(structure, pseudoknots):
    """Creates domain-strings for np-scripting language. Adds
    UG + CG pairs around pseudo-knots, and adds AA + A padding.

    Args:
        structure -- str: The dotbracket secondary structure
        pseudoknots -- list<str>: A list of sequences of the pseudoknots
    Returns:
        str -- A domain following the np-script notation
    """
    t = structure
    for p in pseudoknots:
        #print(p)
        t = re.sub("\(\(\.\.(\[|\])+.\)\)", "UGAA" + p + "ACG", t, count=1)
    t = re.subn("[.()]", "N", t)[0]
    return t


def read_primary(file):
    """Reads the primary structure from NUPACK output

    Args:
        file -- str: The path to the NUPACK output file
    Returns:
        str -- The primary structure
    """
    pattern = ".*Sequence : .*"
    matcher = re.compile(pattern)
    t = open(file)
    matches = [x for x in t if "sequence:" in x]
    sequences = [re.sub(".*sequence: ", "", s)[:-1] for s in matches]
    #TODO: create multiple strands
    return sequences[0].replace(" ", "")


def generate_np(secondary_structure, pseudoknots, GC_content = 1.0, ignore_conflicts = False, ignore_primary = False, restriction=""):
    """Generates the contents of a temporary np file for NUPACK

    args:
        secondary_structure --str: The dotbracket secondary structure sans anglebrackets
        pseudoknots -- list<str>: list of the pseudoknot sequences

    KWargs:
        GC_content -- float: the maximum proportion of G's and C's in the primary structure
        ignore_conflicts -- bool: ignore the conflicts between restrictions and given pseudoknots
        ignore_primary -- bool: ignore the given primary structure as a restriction
        restriction -- str: a list of restrictions for the generator

    returns:
        str -- np contents
        str -- the primary structure in np notation
    """
    simplified = secondary_structure.replace("[", ".").replace("]", ".").replace("A", ".")
    domain = []
    t_d = create_domain(secondary_structure, pseudoknots)
    #print(restriction)
    #print(t_d)
    for i, c in enumerate(restriction):
        #print(i, c)
        if  not ignore_primary and c != "N":
            domain.append(c)
        else:
            domain.append(t_d[i])
    domain = "".join(domain)
    print(domain)
    #print(structure)
    #print(domain)
    print(simplified)
    contents = []
    for l in open(NP_TEMPLATE):
        key, val = l.split("=")
        if "domain" in key:
            val = domain
        elif "structure" in key:
            val = simplified
        if ignore_conflicts and "prevent" in key:
            continue
        contents.append("{}={}".format(key, val))
    #print("\n".join(contents))
    contents = "\n".join(contents)
    return contents, domain


def split_to_domains(primary_structure, secondary_structure):
    """ Splits the primary structure into domains.

    Args:
        primary_structure -- str: IUPAC primary structure
        secondary_structure -- str: dotbracket notation of the secondary structure
    Returns:
        list<str>
    """
    pairs = {}
    stack = []
    p1 = "..[[[[[[."
    p2 = "..]]]]]]."
    ss = secondary_structure.replace(p1, "P1").replace(p2, "P2")
    #print(ss)

    for i, sym in enumerate(ss):
        if sym == "(":
            stack.append(i)
        elif sym == ")":
            pairs[i] = stack.pop()
            pairs[pairs[i]] = i

    domains = []
    for i, sym in enumerate(ss):
        if sym == "(":
            if pairs.get(i - 1) == None:
                domains.append([sym])
            elif pairs[i] == pairs[i - 1] - 1:
                domains[-1].append(sym)
        elif sym == ")":
            if pairs.get(i - 1) == None or pairs[i] != pairs[i - 1] - 1:
                domains.append([sym])
            else:
                domains[-1].append(sym)
        else:
            if len(domains) < 1 or domains[-1][0] in "()":
                domains.append([sym])
            else:
                domains[-1].append(sym)
    for i, d in enumerate(domains):
        if "".join(d) == "P1":
            domains[i] = p1
        elif "".join(d) == "P2":
            domains[i] = p2
    #print(domains)

    p_domains = []
    i = 0
    for t in domains:
        p_domains.append(primary_structure[i : i + len(t)])
        i += len(t)

    return p_domains

if __name__ == "__main__":
    main(sys.argv[1:])
