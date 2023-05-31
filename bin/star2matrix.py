#!/usr/bin/env python3

"""
Loads STAR ReadsPerGene output from multiple samples, and returns sample
metadata and gene count matrices. Expects STAR input files for each
sample in separate folders.
"""

import os
import sys
import numpy as np
import pandas as pd

files = list()

# Set paramters and defaults based on flags
query = ""
stranded = False
combine = False

for i in range(0, len(sys.argv)):
    # STAR inputs
    if sys.argv[i] == "-i":
        files = sys.argv[i + 1 :]

    # Output
    if sys.argv[i] == "-o":
        out = sys.argv[i + 1]

    # Filter options
    if sys.argv[i] == "-q":
        query = sys.argv[i + 1]
    if sys.argv[i] == "-w":
        filter = 0
    if sys.argv[i] == "-b":
        filter = 1

    # Input types
    if sys.argv[i] == "-p":
        stranded = True

    # Combine
    if sys.argv[i] == "-C":
        combine = True

    # Get target
    if sys.argv[i] == "-f":
        strand = 0
    if sys.argv[i] == "-r":
        strand = 1


# Parameter debugging
print(f"Stranded: {stranded}\nCombine: {combine}")


# Parameter logic checks
if (stranded is False) and (combine is True):
    print(
        "Cannot combine counts from non-pair-ended reads, please enable \
        pair-ended data with -p, or disable column combination with -C"
    )
    exit(1)


if len(files) == 0:
    print("Could not find -i flag, or a list of files following the flag")
    exit(1)

# print(f"{files}")


# Set column for strands
if (stranded is True) and (combine is False):
    if strand == 0:  # Set for forward strand
        col = 2
    elif strand == 1:  # Set for reverse strand
        col = 3
elif (stranded is True) and (combine is True):
    col = (2, 3)
else:
    print("Cannot process STAR non-stranded output, not-yet-implemented")
    sys.exit(1)


# Create empty matrix
with open(files[0], "r") as f:
    x = len(f.readlines()) - 4
y = len(files)
matrix = np.zeros((x, y))

print(f"Created empty matrix with dimensions ({x}, {y})")


# Add sample data to each matrix
ids = list()
genes = list()

for i in range(0, len(files)):
    file = files[i]
    id = os.path.dirname(file).split("/")[-1]

    if len(query) != 0:
        if filter == 0:
            id = id.replace(id, query)
        elif filter == 1:
            id = id.replace(query, "")
        else:
            print("Please specify a filter with -w or -b (whitelist/blacklist)")
            exit(1)

    # print(f"ID: {id}")
    ids.append(id)

    with open(file, "r") as f:
        buffer = f.readlines()

        # Remove header
        for j in range(0, 4):
            buffer.pop(0)

        # Add gene count to each matrix
        for j in range(0, len(buffer)):
            line = buffer[j]
            x = line.split("\t")
            gene = x[0]

            if i == len(files) - 1:
                genes.append(gene)

            if stranded is True:
                if combine is True:
                    count = int(x[col[0]]) + int(x[col[1]])
                else:
                    count = int(x[col])
            elif stranded is False:
                count = int(x[col])

            # print(f"Gene: {gene}\tCount: {count}")
            matrix[j, i] = count

        f.close()


# Output matrix
print(f"Number of samples: {len(ids)}")
print(f"Number of genes: {len(genes)}")

df = pd.DataFrame(data=matrix, index=genes, columns=ids)
df[ids] = df[ids].astype(int)
pd.DataFrame.to_csv(df, out, sep="\t")


# Fix header
with open(out, "r") as f:
    buffer = f.readlines()
    f.close()

with open(out, "w") as f:
    buffer[0] = buffer[0][1:]
    for i in buffer:
        f.write(f"{i.rstrip()}\n")
    f.close()
