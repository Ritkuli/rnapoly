#!/usr/bin/python3
# -*- coding: utf-8 -*-
""" Converts a SNAC file to OXDNA input files

Example:
    Show help:
        $ python snac2ox.py -h
    Run with ..:

Author: Antti Elonen
"""
from tempfile import NamedTemporaryFile
from lib.utils import read_snac, remove_file
from shutil import copyfile
import re, os, argparse, sys, math, subprocess, math

OXDNA_EXE = "./lib/oxdna/oxDNA"
RELAXATION_INPUT = "resources/templates/relaxation_input"
GENERIC_INPUT = "resources/templates/generic_input"
BASE_STEPS = 1500
TRAP = {
    "type": "mutual_trap",
    "particle": "n",
    "ref_particle": "n",
    "stiff": "100",
    "r0": "1.2"
}
verbose = True


def main(argv):
    args = interpret_arguments(argv)

    d, path = read_snac(args.input)
    primary_structure = d["primary_structure"]
    positions = d["positions"]
    secondary_structure = d["secondary_structure"]
    pseudoknot_numbering = d["pseudoknot_numbering"]

    bases = get_bases(positions, secondary_structure, pseudoknot_numbering)
    top = generate_top(primary_structure, args.output)
    conf = generate_conf(bases, args.output)
    generated_files = [top, conf]

    input_parameters = read_template(GENERIC_INPUT)
    input_parameters["topology"] = top
    input_parameters["conf_file"] = conf

    if not args.no_relax:
        if args.no_traps:
            external = None
        else:
            external = generate_external(secondary_structure, pseudoknot_numbering, only_pseudoknots = args.only_pseudoknots)
        relaxed_conf = relax_structure(conf, top, external = external)
        if args.replace:
            copyfile(relaxed_conf, conf)
        else:
            rx = conf + ".relaxed"
            copyfile(relaxed_conf, rx)
            generated_files.append(rx)
            input_parameters["conf_file"] = rx
        if not args.no_clean:
            remove_file(relaxed_conf)
            remove_file(external)

    if args.generate_input:
        input_file = args.output + ".input"
        with open(input_file, "w") as f:
            content = "\n".join([p + "=" + str(input_parameters[p]) for p in input_parameters])
            f.write(content)
        generated_files.append(input_file)
    print("**************************")
    print("Generated files: {}".format(", ".join(generated_files)))


def interpret_arguments(argv):
    """Interprets commandline arguments.

    Args:
        argv (list<str>): a list of commandline arguments

    Returns
        namespace: an object with the arguments and their values
    """
    parser = argparse.ArgumentParser()
    g = parser.add_argument_group()

    g.add_argument("input", type=str, default = None,
                    help="The name of the shape")
    g.add_argument("-o", "--output", default=None, type=str,
                    help="The output file name.")
    g.add_argument("-x", "--no_relax", action = "store_true",
                    help="Do not run the relaxation procedure.")
    g.add_argument("-r", "--replace", action = "store_true",
                    help="Replace the initial configuration with the relaxed one.")
    g.add_argument("-c", "--no_clean", action = "store_true",
                    help="Do not clean up the temporary files.")
    g.add_argument("-t", "--no_traps", action = "store_true",
                    help="Do not add traps to the relaxation procedure.")
    g.add_argument("-p", "--only_pseudoknots", action = "store_true",
                    help="Add traps only to the pseudoknots.")
    g.add_argument("-g", "--generate_input", action = "store_true",
                    help="Also generate an input file.")

    args = parser.parse_args(argv[1:])
    if not args.output:
        args.output = args.input.split(".")[0]
    return args

def relax_structure(conf, top, external = None, relaxation_input = RELAXATION_INPUT):
    """ Tries to relax the structure with OXDNA. Generates a new configuration file.

    Args:
        conf -- configuration file
        top -- topology file

    KWargs:
        external -- an external forces file
        relaxation_input -- input parameters template

    Returns:
        str -- path to the generated conf
    """
    min_dt = 0.000000003
    max_dt = 0.003
    dt_multiplier = 5
    input_parameters = read_template(relaxation_input)

    last_conf_name = input_parameters["lastconf_file"]
    copyfile(conf, last_conf_name)
    input_parameters["external_forces_file"] = external
    input_parameters["topology"] = top
    input_parameters["conf_file"] = last_conf_name
    input_parameters["steps"] = int(input_parameters["steps"])
    if external != None: input_parameters["external_forces"] = 1
    print("Starting relaxation procedure.")
    for relax_strength in [10, 100]:
        i = 0
        input_parameters["relax_strength"] = relax_strength
        print("-- Starting a cycle with relax strength {}".format(relax_strength))
        while min_dt * dt_multiplier ** i <= max_dt:
            iterations = int(math.log(max_dt / min_dt, dt_multiplier))
            input_parameters["dt"] = min_dt * dt_multiplier ** i
            if (run_relaxation(input_parameters) == 0):
                #copyfile(last_conf_name, last_conf_name + ".bak")
                i += 1
                print(".. Step {} of {} complete.".format(i, iterations + 1))
            else:
                print(".. Failed step {}. Reverting back to start.".format(i + 1))
                copyfile(conf, last_conf_name)
                if i > 0:
                    input_parameters["steps"] *= 2
                    print(".. Increasing step count to {}".format(input_parameters["steps"]))
                else:
                    min_dt /= dt_multiplier
                    print(".. Increasing step count by one.")
                i = 0
    return last_conf_name

def generate_external(secondary_structure, pseudoknot_numbering, only_pseudoknots = False, five_to_three = False):
    """Generates an external forces file for OXDNA.

    Args:
        secondary_structure
        pseudoknot_numbering

    KWArgs:
        only_pseudoknots -- Create external forces only for pseudoknots
        five_to_three -- read from five prime to three prime

    Returns:
        str -- path to the generated file
    """
    path = ".external.conf.temp"
    if five_to_three:
        ss = secondary_structure
        pn = pseudoknot_numbering.split(" ")
    else:
        ss = secondary_structure[::-1]
        pn = pseudoknot_numbering[::-1].split(" ")
    n = None
    d = {}
    closing = False
    stack = []
    f = open(path, "w")
    for i, c in enumerate(secondary_structure):
        if c in "P":
            if n == None:
                n = pn.pop(0)
                if n in d: closing = True
            if closing:
                f.writelines(get_traps(d[n].pop(), i))
            else:
                d.setdefault(n, []).append(i)
        else:
            n = None
            closing = False
            if not only_pseudoknots:
                if c == "(":
                    stack.append(i)
                elif c == ")":
                    f.writelines(get_traps(stack.pop(), i))
    f.close()
    return path

def get_traps(i1, i2):
    """ Generates kinetic traps for the given indices of bases

    args:
        i1 -- first index
        i2 -- second index

    returns:
        str -- the trap
    """
    result = []
    for i in [i1, i2]:
        result.append("{")
        for k in TRAP:
            if k == "particle":
                key, val = k, i
            elif k == "ref_particle":
                key, val = k, i1 + i2 - i
            else:
                key, val = k, TRAP[k]
            result.append("{} = {}".format(key, val))
        result.append("}\n")
    return "\n".join(result)


def read_template(template):
    """ Reads an OXDNA template into a dictionary

    args:
        template -- path of the template

    returns:
        dictionary -- argument => value
    """
    content = {}
    with open(template) as f:
        for l in f:
            line = l.strip()
            if line.startswith("#") or len(line) < 1:
                continue
            p, k = line.split("=")
            content[p.strip()] = k.strip()
    return content

def run_relaxation(input_params):
    """ Runs relaxation for the input parameters dictionary.

    returns:
        bool -- 0 if successful
    """
    t = NamedTemporaryFile(suffix = ".dat", mode = "w")
    content = "\n".join([p + "=" + str(input_params[p]) for p in input_params])
    t.write(content)
    t.seek(0)
    x = ""
    cmd = [OXDNA_EXE, t.name]
    kwargs = {"stderr": subprocess.PIPE, "stdout": subprocess.PIPE}
    #if verbose: kwargs["stdout"] = None
    with subprocess.Popen(cmd, **kwargs) as sp:
        while sp.poll() == None:
            try:
                a, b = sp.communicate(timeout=1)
            except subprocess.TimeoutExpired as e:
                if "displacement" in str(sp.stderr.readline()):
                    sp.kill()
                    return 1
    return 0

def get_direction(pos1, pos2):
    """ Gets the normalized direction vector from pos1 to pos2
    Args:
        pos1
        pos2

    Returns:
        direction
    """
    d = [t[1] - t[0] for t in zip(pos1, pos2)]
    return normalize(d)

def normalize(v):
    """ Normalizes v
    Args:
        v -- vector

    Returns:
        vector
    """
    l = sum([t ** 2 for t in v]) ** -0.5
    return [l * t for t in v]

def cross(v1, v2):
    """ Calculates the cross product of v1 and v2
    Args:
        v1
        v2

    Returns:
        v1 x v2
    """
    n = [v1[1] * v2[2] - v1[2] * v2[1],
         v1[2] * v2[0] - v1[0] * v2[2],
         v1[0] * v2[1] - v1[1] * v2[0]]
    return normalize(n)

def distance(v1, v2):
    """ Calculates the euclidian distance between v1 and v2
    Args:
        v1
        v2

    Returns:
        float
    """
    n = sum([(x - y)**2 for (x, y) in zip(v1, v2)])
    return n ** 0.5


def generate_conf(base_list, output, five_to_three = False, size = 100):
    """Generates a confiugration file for oxrna

    Args:
        base_list -- bases in 5'-3'-order
        output -- output path

    KWargs:
        five_to_three -- read from five prime to three prime
        size -- the size of the bounding cube edge

    Returns:
        str -- path of the generated file
    """
    path = output.split(".top")[0]
    if not path.endswith(".conf"):
        path += ".conf"
    print("Generating the initial confiugration file {}".format(path))

    bases = base_list[:] if five_to_three else base_list[::-1]

    T = 0
    L = " ".join([str(t) for t in 3 * [size]])
    E = "0 0 0"

    f = open(path, "w")
    f.write("t = {}\nb = {}\nE = {}\n".format(T, L, E))
    direction = [1, 0, 0]
    base_versor = [0, 1, 0]
    last = None
    for i in range(len(bases)):
        base = bases[i]
        if base.pair:
            base_versor = [(1.0 / 1.0)*x for x in get_direction(base.position, base.pair.position)]
        else:
            base_versor = [x for x in base_versor]
        base_position = [((1.0 / 1.0)*x + 0 * y) for x, y in zip(base.position, base_versor)]
        if last:
            direction = get_direction(last, base.position)
            print(distance(last, base_position))
        last = base_position
        norm = [(1.0 / 1.0) * x for x in cross(base_versor, direction)]

        r = " ".join([str(t) for t in base_position])
        b = " ".join([str(t) for t in base_versor])
        n = " ".join([str(t) for t in norm])
        v = "0 0 0"
        l = "0 0 0"

        f.write("{} {} {} {} {}\n".format(r, b, n, v, l))
    #sys.exit(0)
    return path

def generate_top(primary_structure, output, five_to_three = False):
    """Generates a topology file for oxrna

    Args:
        primary_structure
        output -- path

    KWArgs:
        five_to_three -- read from 5' to 3'

    Returns:
        str -- output path
    """
    path = output.split(".conf")[0]
    if not path.endswith(".top"):
         path += ".top"
    print("Generating a topology file {}".format(path))

    ps = primary_structure if five_to_three else primary_structure[::-1]
    f = open(path, "w")
    strands = 1
    for s in range(1, strands + 1):
        num = len(ps)
        f.write("{} {}\n".format(num, strands))
        for i, b in enumerate(ps):
            prev = i - 1
            next = i + 1 if i < num - 1 else -1
            f.write("{} {} {} {}\n".format(s, b, prev, next))
    f.close()
    return path


def db_to_dict(s_str, i = 0, d = {}):
    """ Converts a dotbracket string to a dictionary

    Args:
        s_str -- secondary_structure

    KWargs:
        i -- start index
        d -- the dictionary so far

    Returns:
        dictionary
    """
    j = i
    while j < len(s_str):
        c = s_str[j]
        if c is "(":
            d = db_to_dict(s_str, j + 1, d)
            j = d[j]
        elif c is ")":
            d[i - 1] = j
            d[j] = i - 1
            if(i != 0): return d # Don't return from the first iteration yet
        else:
            d[j] = None
        j = j + 1
    return d


class Base():
    """ A representation of a base containing its position, pseudoknot number
    and its pair
    """

    def __init__(self, position, p_num = None):
        self.position = [float(t) for t in position]
        self.p_num = p_num
        self.pair = None

    def set_pair(self, other):
        """ Pairs this pair with the other
        """
        self.pair = other
        other.pair = self


def get_bases(positions, secondary_structure, pseudoknot_numbering):
    """Generates a list of bases

    Args:
        positions
        secondary_structure
        pseudoknot_numbering

    Returns:
        list
    """
    d = db_to_dict(secondary_structure)
    pos_l = positions.strip("[ ]").split(",")
    p_num = None
    bases = []
    pn = pseudoknot_numbering.split(" ")
    for i in range(len(secondary_structure)):
        if secondary_structure[i] == "P":
            if not p_num:
                p_num = pn.pop(0)
        else:
            p_num = None
        pos = pos_l[i].strip().split(" ")
        b = Base(pos, p_num)
        other_i = d[i]
        if other_i and other_i < i:
            b.set_pair(bases[other_i])
        bases.append(b)
    for b in bases:
        #print(b.position)
        pass
    return bases

if __name__ == "__main__":
    main(sys.argv)
