#!/usr/bin/env python3

"""
Prints a given value from a TOML key to stdout, given the path to the TOML file
and a dot (.) delimited list of keys.
"""

import sys
try:
    import toml
except ImportError:
    exit(1)

file = sys.argv[1]
keys = sys.argv[2].split(".")

with open(file, "r") as f:
    tf = toml.load(file)
    f.close()

for i in range(0, len(keys)):
    key = keys[i]
    if i == 0:
        last_key = tf[key]
    else:
        last_key = last_key[key]
print(last_key)
