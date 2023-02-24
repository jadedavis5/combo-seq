#!/usr/bin/env python3

"""
Generate Nextflow configurations from TOML files.

Takes
"""
import sys
try:
    import toml
except ImportError:
    exit(1)

args = sys.argv

def load_toml(file):
    with open(file, "r") as f:
        return toml.load(f)
        print("Closed TOML file")
        f.close()

def nf_conf_str(x, conf):
    string = "\n"
    for i in conf[x]:
        val = conf[x][i]
        if isinstance(val, str):
            string += f"\t{i} = \"{val}\"\n"
        else:
            string += f"\t{i} = {val}\n"
    return f"{x} {{{string}}}\n"

# Load vars from workflow config
settings = load_toml("settings.toml")
params = settings["files"]["params"]
nf_conf = settings["files"]["config"]

toml_common = settings["files"]["common"]
toml_in = args[1]

conf = load_toml(toml_in)
conf_common = load_toml(toml_common)

# Write yaml file
with open(params, "w") as f:
    for i in conf_common["data"]:
        f.write(f"{i}: \"{conf_common['data'][i]}\"\n")
    f.close()

# Write Nextflow configuration
with open(nf_conf, "w") as f:
    f.write(nf_conf_str("process", conf))
    f.close()
