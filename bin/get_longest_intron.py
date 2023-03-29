#!/usr/bin/env python3

"""
Get length of introns from a GFF file containing intron annotations.
"""

import os
import sys
from pprint import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer as GE

for i, arg in enumerate(sys.argv):
    if arg == "-a":
        infile = sys.argv[i + 1]

introns = set()
with open(infile, "r") as f:
    # Shows annotation levels
    # pprint(GE().parent_child_map(f))

    # Show high-level summary
    # pprint(GE().available_limits(f))

    # Get intronic info
    for i in GFF.parse(f, limit_info = dict(gff_type = ["intron"])):
        introns.add((i.features[0].location.start.position, i.features[0].location.end.position))
    f.close()

ranges = set()
for i in introns:
    ranges.add(i[1] - i[0])

print(max(ranges))
