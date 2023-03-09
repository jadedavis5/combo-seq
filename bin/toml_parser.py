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
    x = toml.load(file)
    f.close()

key = x[keys[0]]
for i in range(1, len(keys)):
    key = key[keys[i]]
print(key)
