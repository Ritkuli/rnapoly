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
from shutil import copyfile, which
import re, os, argparse, sys, math, subprocess, math, tempfile, numpy as np

script_path = os.path.abspath(os.path.dirname(sys.argv[0]))
work_path = os.getcwd()

OXDNA_EXE = which("oxDNA")
if not OXDNA_EXE: OXDNA_EXE = script_path + "/lib/oxdna/oxDNA"
RELAXATION_INPUT = script_path + "/resources/templates/relaxation_input"
BURN_IN_INPUT = script_path + "/resources/templates/burn_in_input"
GENERIC_CPU_INPUT = script_path + "/resources/templates/generic_input"
GENERIC_GPU_INPUT = script_path + "/resources/templates/oxrna_cuda_input"
SEQ_DEP_FILE = script_path + "/resources/rna_sequence_dependent_parameters.txt"
BASE_STEPS = 1000
TRAP = {
    "type": "mutual_trap",
    "particle": "n",
    "ref_particle": "n",
    "stiff": "15",
    "r0": "1.2"
}
verbose = False


"""
pairs: (1, 3); (2, 4)
5' -> 3':1 -> 2 -> 4 -> 3
"""
magic_value = 1.15
# Backbone coordinates of the reference RNA helix:
reference_helix_bb = np.matrix(([[0.34790706782241565, 1.425120976291932, 1.9743945943644676, 1], [-0.10455815001983909, 0.9425201267314649, 2.3030945943644676, 1], [-0.22459207165999268, 0.2919663541582757, 2.6317945943644676, 1], [0.025853467663527308, -0.32032882399232404, 2.9604945943644676, 1]])).transpose()
# Inverse of the above:
reference_helix_bb_inv = np.linalg.inv(reference_helix_bb)
# OXDNA input parameters of the reference helix:
reference_helix_input_1 = np.matrix(([[0.6983741093562843, 1.264664636888758, 2.0812899447957705, 1.0], [0.8761676038346716, -0.4011408485079351, 0.2672383760782569, 1.0], [-0.24298278125345144, 0.11124620291623981, 0.963630453208623, 1.0]])).transpose()
reference_helix_input_2 = np.matrix(([[0.27704862836922695, 0.9968308139273313, 2.4099899447957704, 1.0], [0.9540169459726651, 0.13577671798966606, 0.26723837607825685, 1.0], [-0.2645723145672288, -0.0376542164104254, 0.963630453208623, 1.0]])).transpose()
reference_helix_input_3 = np.matrix(([[0.0671933237617629, 0.5438287512618107, 2.7386899447957704, 1.0], [0.7294634885543889, 0.6296559927588375, 0.2672383760782569, 1.0], [-0.20229812937164504, -0.17461906110638953, 0.963630453208623, 1.0]])).transpose()
reference_helix_input_4 = np.matrix(([[0.13532780179762505, 0.04925033427023945, 3.0673899447957704, 1.0], [0.2736858353352443, 0.9239478956564087, 0.2672383760782569, 1.0], [-0.07589979949998266, -0.2562334288979607, 0.963630453208623, 1.0]])).transpose()


def main(argv):
    args = interpret_arguments(argv)

    d, path = read_snac(args.input)
    primary_structure = d["primary_structure"]
    positions = d["positions"]
    secondary_structure = d["secondary_structure"]
    pseudoknot_numbering = d["pseudoknot_numbering"]

    five_to_three = False

    bases = get_bases(positions, primary_structure, secondary_structure, pseudoknot_numbering, five_to_three)

    top = os.path.join(work_path, generate_top(bases, args.output))
    conf = os.path.join(work_path, generate_conf(bases, args.output))
    generated_files = [top, conf]

    input_parameters = {}
    input_parameters["topology"] = top
    input_parameters["conf_file"] = conf
    input_parameters["seq_dep_file"] = SEQ_DEP_FILE

    burn_in = None
    external = None
    relaxed_conf = None

    with tempfile.TemporaryDirectory() as temp_dir_path:
        os.chdir(temp_dir_path)
        if not args.no_traps:
            external = generate_external(bases, only_pseudoknots = args.only_pseudoknots)

        if not args.no_relax:
            relaxed_conf = relax_structure(conf, top, external = external, steps = args.relax_steps, iterations = args.relax_iterations)
            if args.replace:
                copyfile(relaxed_conf, conf)
            else:
                rx = conf + ".relaxed"
                copyfile(relaxed_conf, rx)
                generated_files.append(rx)
                input_parameters["conf_file"] = rx

        if not args.no_burn_in:
            burn_in = burn_in_structure(input_parameters["conf_file"], top, external = external, steps = args.burnin_steps)
            if args.replace:
                copyfile(burn_in, conf)
            else:
                bi = conf + ".burnin"
                copyfile(burn_in, bi)
                generated_files.append(bi)
                input_parameters["conf_file"] = bi

        if not args.absolute_paths:
            out_dir = os.path.dirname(input_parameters["conf_file"])
            input_parameters["conf_file"] = os.path.basename(input_parameters["conf_file"])
            input_parameters["topology"] = os.path.basename(input_parameters["topology"])
            input_parameters["seq_dep_file"] = os.path.relpath(input_parameters["seq_dep_file"], out_dir)
        os.chdir(work_path)

    if args.generate_cpu_input:
        input_parameters = {**read_template(GENERIC_CPU_INPUT), **input_parameters}
        input_file = generate_input(input_parameters, args.output + "_cpu")
        generated_files.append(input_file)

    if args.generate_gpu_input:
        input_parameters = {**read_template(GENERIC_GPU_INPUT), **input_parameters}
        input_file = generate_input(input_parameters, args.output + "_gpu")
        generated_files.append(input_file)


    print("**************************")
    print("Generated files: ")
    # Make output human readable
    for f in generated_files:
        f = os.path.relpath(f, work_path)
        print(f)


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
    g.add_argument("-b", "--no_burn_in", action = "store_true",
                    help="Do not run the burn-in simulation.")
    g.add_argument("-xs", "--relax_steps", default=500, type=int,
                    help="The number of steps in the relaxation process. Default = 500")
    g.add_argument("-xi", "--relax_iterations", default=7, type=int,
                    help="The number of iterations in the relaxation process. Default = 7")
    g.add_argument("-bs", "--burnin_steps", default=10000, type=int,
                    help="The number of steps in the burn-in process. Default = 10000")
    g.add_argument("-r", "--replace", action = "store_true",
                    help="Replace the initial configuration with the relaxed one.")
    g.add_argument("-t", "--no_traps", action = "store_true",
                    help="Do not add traps to the relaxation procedure.")
    g.add_argument("-p", "--only_pseudoknots", action = "store_true",
                    help="Add traps only to the pseudoknots.")
    g.add_argument("-g", "--generate_cpu_input", action = "store_true",
                    help="Also generate an input file for CPU simulations.")
    g.add_argument("-G", "--generate_gpu_input", action = "store_true",
                    help="Also generate an input file for GPU simulations.")
    g.add_argument("-a", "--absolute_paths", action = "store_true",
                    help="Use absolute instead paths of relative for the input-file.")

    args = parser.parse_args(argv[1:])
    if not args.output:
        args.output = "".join(args.input.split(".")[0:-1])

    return args

def relax_structure(conf, top, external = None, relaxation_input = RELAXATION_INPUT, steps = 500, iterations = 6):
    """ Tries to relax the structure with OXDNA. Generates a new configuration file.

    Args:
        conf -- str: configuration file
        top -- str: topology file

    KWargs:
        external -- str: an external forces file
        relaxation_input -- str: input parameters template

    Returns:
        str -- path to the generated conf
    """
    dt_multiplier = 10
    min_dt = 0.003 * dt_multiplier ** (-iterations + 1)
    max_dt = 0.003
    input_parameters = read_template(relaxation_input)

    last_conf_name = input_parameters["lastconf_file"]
    copyfile(conf, last_conf_name)
    input_parameters["external_forces_file"] = external
    input_parameters["topology"] = top
    input_parameters["conf_file"] = last_conf_name
    input_parameters["steps"] = steps
    if external != None: input_parameters["external_forces"] = 1
    print("Starting relaxation procedure.")
    for relax_strength in [1, 10, 100]:
        i = 0
        input_parameters["relax_strength"] = relax_strength
        print("-- Starting a cycle with relax strength {}".format(relax_strength))
        while min_dt * dt_multiplier ** i <= max_dt:
            iterations = math.ceil(math.log(max_dt / min_dt, dt_multiplier))
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

def burn_in_structure(conf, top, external = None, burn_in_input = BURN_IN_INPUT, steps = 500):
    """ Burns in the RNA strand by forcing basepairs

    Args:
        conf -- str: path to the configuration file
        top -- str: path to the topology file
    Kwargs:
        external -- str: path ot the external forces file
        burn_in_input -- str: path to the burn-in-input parameters file
    Returns:
        str -- path to the created file
    """
    input_parameters = read_template(burn_in_input)
    last_conf_name = input_parameters["lastconf_file"]
    copyfile(conf, last_conf_name)
    input_parameters["external_forces_file"] = external
    input_parameters["seq_dep_file"] = SEQ_DEP_FILE
    input_parameters["topology"] = top
    input_parameters["conf_file"] = last_conf_name
    input_parameters["steps"] = steps
    if external != None: input_parameters["external_forces"] = 1
    print("Starting burn-in procedure.")
    if run_relaxation(input_parameters) != 0:
        print("Error! Try relaxing the structure or increasing relaxation steps / iterations.")
    return last_conf_name

def generate_external(bases, only_pseudoknots = False):
    """Generates an external forces file for OXDNA.

    Args:
        bases -- list<Base>
    KWArgs:
        only_pseudoknots -- bool, create external forces only for pseudoknots
    Returns:
        str -- path to the generated file
    """
    path = ".external.conf.temp"

    f = open(path, "w")
    for i, b in enumerate(bases):
        if b.p_num:
            f.writelines(get_traps(bases.index(b.pair), i))
        else:
            if not only_pseudoknots and b.pair:
                f.writelines(get_traps(bases.index(b.pair), i))
    f.close()
    return path

def generate_input(input_parameters, output):
    """Generates an OXDNA-simulation input file.

    Args:
        input_parameters <dict<param, val>>
        output -- str: output path

    Returns:
        str -- path to the generated file
    """
    input_file = os.path.join(work_path, output + ".input")

    with open(input_file, "w") as f:
        content = "\n".join([p + "=" + str(input_parameters[p]) for p in input_parameters])
        f.write(content)

    return input_file

def get_traps(i1, i2):
    """ Generates kinetic traps for the given indices of bases

    args:
        i1 -- int, first index
        i2 -- int, second index

    returns:
        str -- the trap in oxdna input format
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
        template -- str: path of the template

    returns:
        dict<arg, val>
    """
    content = {}
    with open(template) as f:
        for l in f:
            line = l.strip() #.replace("$rpoly", script_path)
            if line.startswith("#") or len(line) < 1:
                continue
            p, k = line.split("=")
            content[p.strip()] = k.strip()
    return content

def run_relaxation(input_params):
    """ Runs relaxation for the input parameters dictionary.

    Args:
        input_params -- dict<param, val>
    Returns:
        bool -- 0 if successful
    """
    t = NamedTemporaryFile(suffix = ".dat", mode = "w")
    content = "\n".join([p + "=" + str(input_params[p]) for p in input_params])
    t.write(content)
    t.seek(0)
    x = ""
    cmd = [OXDNA_EXE, t.name]
    kwargs = {"stderr": subprocess.PIPE, "stdout": subprocess.PIPE}
    if verbose: kwargs = {}
    with subprocess.Popen(cmd, **kwargs) as sp:
        while sp.poll() == None:
            try:
                a, b = sp.communicate(timeout=1)
            except subprocess.TimeoutExpired as e:
                print("...")
        if sp.poll() == 1:
            return 1
    return 0

def get_direction(pos1, pos2):
    """ Gets the normalized direction vector from pos1 to pos2
    Args:
        pos1 -- Vector
        pos2 -- Vector

    Returns:
        direction -- Vector
    """
    d = [t[1] - t[0] for t in zip(pos1, pos2)]
    return normalize(d)

def normalize(v):
    """ Normalizes v
    Args:
        v -- Vector

    Returns:
        Vector
    """
    l = sum([t ** 2 for t in v]) ** -0.5
    return [l * t for t in v]

def cross(v1, v2):
    """ Calculates the cross product of v1 and v2
    Args:
        v1 -- Vector
        v2 -- Vector

    Returns:
        Vector -- v1 x v2
    """
    n = [v1[1] * v2[2] - v1[2] * v2[1],
         v1[2] * v2[0] - v1[0] * v2[2],
         v1[0] * v2[1] - v1[1] * v2[0]]
    return normalize(n)

def distance(v1, v2):
    """ Calculates the euclidian distance between v1 and v2
    Args:
        v1 -- Vector
        v2 -- Vector

    Returns:
        float
    """
    n = sum([(x - y)**2 for (x, y) in zip(v1, v2)])
    return n ** 0.5


def generate_conf(bases, output, size = 100):
    """Generates a confiugration file for oxrna

    Args:
        bases -- list<Base>: bases in 5'-3'-order
        output -- str: output path
    KWargs:
        size -- float: the size of the bounding cube edge
    Returns:
        str -- path of the generated file
    """
    path = output.split(".top")[0]
    if not path.endswith(".conf"):
        path += ".conf"
    print("Generating the initial confiugration file {}".format(path))

    T = 0
    L = " ".join([str(t) for t in 3 * [size]])
    E = "0 0 0"

    f = open(path, "w")
    f.write("t = {}\nb = {}\nE = {}\n".format(T, L, E))

    i = 0
    stack = []
    last_input = None
    transform = None # used if next base is unpaired, i.e., at the end of helix
    for i, base in enumerate(bases):
        if not base.pair:
            stack.insert(0, base)
            continue
        try:
            if all([bases[i + j].pair for j in range(4)]):
                n_position = [bases[i + j].position for j in range(4)]
                helix_bb = np.matrix((*n_position,)).transpose()
                #print("----")
                #print(helix_bb)
                transform = helix_bb * reference_helix_bb_inv
        except IndexError:
            pass
        if i + 3 < len(bases) and all([bases[i + j].pair for j in range(4)]):
            input_mats = transform, reference_helix_input_1
        elif i + 2 < len(bases) and all([bases[i + j].pair for j in range(3)]):
            input_mats = transform, reference_helix_input_2
        elif i + 1 < len(bases) and all([bases[i + j].pair for j in range(2)]):
            input_mats = transform, reference_helix_input_3
        else:
            input_mats = transform, reference_helix_input_4
            #base.pair.input_mats = transform, reference_helix_input_4
        #print(np.linalg.cond(transform))
        #print("===")
        input_pos = (input_mats[0] * input_mats[1])[:,0]
        #print(input_pos)
        input_dir_1 = (np.linalg.inv(input_mats[0]).transpose() * input_mats[1])[:,1]
        input_dir_2 = (np.linalg.inv(input_mats[0]).transpose() * input_mats[1])[:,2]
        #print(np.linalg.norm(input_dir_1))
        #print(":::", input_mats[0] * input_mats[1])
        current_input = (input_pos, input_dir_1, input_dir_2)
        #print(current_input)
        #print("---")
        #process stack
        f.write(process_unpaired(stack, last_input, current_input))
        stack.clear()

        last_input = current_input
        print(i, base)
        f.write(extend_conf(*current_input))
    if len(stack) > 0:
        f.write(process_unpaired(stack, last_input, current_input))
    #sys.exit(0)
    print(len(bases))
    f.close()
    return path

def process_unpaired(stack, input1, input2):
    """ Handles the transformation of the unpaired nucleotides

    Args:
        stack -- list<Base>
        input1 -- Vector: last transformed vector
        input2 -- Vector: the next transformed vector
    Returns:
        str -- transformed bases in oxdna input format
    """
    r = []
    size = len(stack) + 1
    for j, b in enumerate(stack):
        if np.all(input1) and not np.all(np.equal(input1, input2)):
            print("Interpolating unpaired base {}.".format(b))
            s_input = [((j + 1) / size) * x + ((size - 1 - j) / size) * y for x, y in zip(input1, input2)]
        else:
            print("Warning: Unable to interpolate unpaired base {}.".format(b))
            s_input = [x.copy() for x in input2]
            s_input[0] = s_input[0] + 0.4 * np.matrix([j, 1, 0, 0]).transpose()
            #print(s_input)
        r.append(extend_conf(*s_input))
    return "".join(r)

def extend_conf(input_pos, input_dir_1, input_dir_2):
    """ Create a new entry for the configuration file

    Args:
        input_pos -- Vector: center of mass
        input_dir_1 -- Vector: base versor
        input dir_2 -- Vector: backbone versor
    Returns:
        str -- the new entry in oxdna input format
    """
    #print(input_dir_1 / input_dir_1[-1])
    r = " ".join([str(t) for t in list((input_pos * magic_value).flat)[:3]])
    b = " ".join([str(t) for t in list((input_dir_1 / np.linalg.norm(input_dir_1[:3]) * 1.0).flat)[:3]])
    n = " ".join([str(t) for t in list((input_dir_2 / np.linalg.norm(input_dir_2[:3]) * 1.0).flat)[:3]])
    v = "0 0 0"
    l = "0 0 0"

    nl = "{} {} {} {} {}\n".format(r, b, n, v, l)
    return nl


def generate_top(bases, output):
    """Generates a topology file for oxrna

    Args:
        bases -- list<Base>
        output -- str: path to the output file
    Returns:
        str -- path to the created file
    """
    path = output.split(".conf")[0]
    if not path.endswith(".top"):
         path += ".top"
    print("Generating a topology file {}".format(path))

    f = open(path, "w")
    strands = 1
    for s in range(1, strands + 1):
        num = len(bases)
        f.write("{} {}\n".format(num, strands))
        for i, b in enumerate(bases):
            prev = i - 1
            next = i + 1 if i < num - 1 else -1
            f.write("{} {} {} {}\n".format(s, b.type, prev, next))
    f.close()
    return path


def db_to_dict(s_str, i = 0, d = {}):
    """ Converts a dotbracket string to a dictionary of indices and their pairs

    Args:
        s_str -- str: secondary_structure in dotbracket notation
    KWargs:
        i -- int: start index
        d -- dict<index1, index2>: the dictionary so far
    Returns:
        dictionary
    """
    j = i
    while j < len(s_str):
        c = s_str[j]
        if c == "(":
            d = db_to_dict(s_str, j + 1, d)
            j = d[j]
        elif c == ")":
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

    def __init__(self, type, position, p_num = None):
        self.type = type
        self.position = np.array(position + [1], dtype=np.float)
        self.p_num = p_num
        self.pair = None

    def set_pair(self, other):
        """ Pairs this pair with the other
        """
        self.pair = other
        other.pair = self

    def dist(self, other):
        return np.linalg.norm(other.position - self.position)

    def translate(self, v):
        self.position += v
        print(self.position)

    def __repr__(self):
        if self.pair:
            return "{} -- {}".format(self.type, self.pair.type)
        else:
            return "{}".format(self.type)


def get_bases(positions, primary_structure, secondary_structure, pseudoknot_numbering, five_to_three = False):
    """Generates a list of bases

    Args:
        positions -- str: serialized 3-tuples of Vectors
        primary_structure -- str: IUPAC primary structure
        secondary_structure -- str: dotbracket secondary structure
        pseudoknot_numbering -- str: list of pseudoknot numbers according to SNAC-format
    Kwargs:
        five_to_three -- bool: read 5'-3'-direction instead of 3'-5'
    Returns:
        list<Base>
    """
    d = db_to_dict(secondary_structure)
    pos_l = positions.strip("[ ]").split(",")
    p_num = None
    bases = []
    pn = pseudoknot_numbering.split(" ")
    for i in range(len(secondary_structure)):
        if secondary_structure[i] in "[]":
            if not p_num:
                p_num = pn.pop(0)
        else:
            p_num = None
        pos = pos_l[i].strip().split(" ")
        b = Base(primary_structure[i], pos, p_num)
        other_i = d[i]
        if other_i != None and other_i < i:
            b.set_pair(bases[other_i])
        bases.append(b)
    #set pseduoknot Pairs
    stacks = {}
    for b in bases:
        if b.p_num:
            stacks.setdefault(b.p_num, []).append(b)
    for i in stacks:
        stack = stacks[i]
        for j in range(int(len(stack) / 2)):
            other = stack[-j - 1]
            stack[j].set_pair(other)

    return bases if five_to_three else bases[::-1]

if __name__ == "__main__":
    main(sys.argv)
