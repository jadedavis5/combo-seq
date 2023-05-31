#!/usr/bin/env bash

# Currently enabled:
# -p for stranded reads
# -C to combine read strands
# -q applies a extact match query string to remove from sample IDs with the -b flag
# -i with the path to two samples

# Using -f or -r without the combine flag (-C) uses only the forward or reverse
# counts column

./bin/star2matrix.py \
    -p -C \
    -q '-wheat-cs' -b \
    -o matrix.tsv \
    -i data/star-wheat/*/ReadsPerGene.out.tab
