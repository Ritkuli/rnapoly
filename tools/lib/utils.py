import threading, subprocess, os
DEVNULL = open(os.devnull, 'wb')


def read_snac(snac):
    """

    Args:
        snac -- The path to the SNAC file

    Returns:
        dictionary -- str => str,  Map of parameters to their values
        path -- str, the absolute path
    """
    if snac.endswith(".snac"):
        path = snac
    else:
        path = snac + ".snac"
    d = {}
    with open(path) as f:
        last = None
        for line in f:
            l = line.split("#")[0]
            if ":" in l:
                key, val = l.split(":")
                last = key.strip()
                val = val.strip()
            else:
                val = l.strip()
            if last:
                t = d.setdefault(last, [])
                t.append(val)
            else:
                print("Unexpected input: ", last)
    for a in d:
        d[a] = (" ".join(d[a])).strip("[ ]")
    return d, path

def write_snac(snac, content, save_over=False):
    """ Writes the contents to a snac file. If not save_over, and if the file
    already exists, only writes the fields that are defined in the input and
    leaves everything else as is.

    Args:
        snac -- The path to the SNAC file
        content -- A dictionary mapping field names to their contents

    KWargs:
        save_over -- ignore previously existing file

    Returns:
        dictionary -- str => str, Map of parameters to their values
    """

    path = snac if snac.endswith(".snac") else snac + ".snac"

    keys = list(content.keys())
    lines = []
    if not save_over and os.path.exists(path):
        with open(path, "r") as f:
            for l in f:
                lines.append(l)
    written = set()
    for i in range(len(lines)):
        l = lines[i]
        if ":" in l:
            key, val = l.split(":")
            key = key.strip()
            val = val.strip()
            if key in content:
                lines[i] = "{}: {}\n".format(key, content[key])
                if "#" in val:
                    lines[i] += "# {}".format(l.split("#")[1])
                written.add(key)
    for key in content:
        if key not in written:
            lines.append("{}: {}\n".format(key, content[key]))

    with open(path, "w") as f:
        f.writelines(lines)

    return path


def remove_file(path):
    """ Removes the file at path

    returns:
        bool -- success
    """
    try:
        os.remove(path)
        return 0
    except FileNotFoundError:
        return 1



class RunCmd(threading.Thread):
    """ Runs a subprocess with the given timeout and command.
    """

    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
        self.p = subprocess.Popen(self.cmd, stdout=DEVNULL, stderr=DEVNULL)
        self.p.wait()

    def run_the_process(self):
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            self.p.kill()
            self.join()
            return -1
        return 1
