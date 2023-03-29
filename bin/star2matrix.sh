#!/usr/bin/env bash

# Currently enabled:
# -p for stranded reads
# -C to combine read strands
# -q with a filter, and -b to remove from the sample IDs
# -i with the path to two samples

./bin/star2matrix.py \
    -p -C \
    -q '-wheat-cs' -b \
    -o matrix.tsv \
    -i data/star/*/ReadsPerGene.out.tab
