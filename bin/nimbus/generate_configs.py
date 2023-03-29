#!/usr/bin/env python3

"""
Generate Nextflow configurations from TOML files.
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
settings = load_toml("conf/nimbus.toml")
params = settings["files"]["params"]
nf_conf = settings["files"]["config"]

toml_common = settings["files"]["common"]
toml_in = args[1]

conf = load_toml(toml_in)
conf_common = load_toml(toml_common)


# Merge dictionaries
x = set()
for k1 in conf:
    v1 = conf[k1]
    # Check for shared keys
    for k2 in conf_common:
        if k2 in conf:
            if k2 == k1:
                v2 = conf_common[k2]
                for i in v2:
                    v1[i] = v2[i]
        else:
            x.add(k2)

for k in x:
    conf[k] = conf_common[k]


# Write YAML file
# with open(params, "w") as f:
#     for i in conf_common:
#         f.write(f"{i}:\n")
#         for j in conf_common[i]:
#             f.write(f"  {j}: \"{conf_common[i][j]}\"\n")
#         f.write("\n")
#     f.close()

with open(params, "w") as f:
    for i in conf:
        f.write(f"{i}:\n")
        for j in conf[i]:
            f.write(f"  {j}: \"{conf[i][j]}\"\n")
        f.write("\n")
    f.close()

with open(params, "a") as f:
    for i in settings:
        f.write(f"{i}:\n")
        for j in settings[i]:
            f.write(f"  {j}: \"{settings[i][j]}\"\n")
        f.write("\n")
    f.close()


# Write Nextflow configuration
with open(nf_conf, "w") as f:
    f.write(nf_conf_str("process", conf))
    f.close()
