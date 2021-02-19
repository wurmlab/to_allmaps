#!/usr/bin/env bash
# Copyright 2019 Anurag Priyam, Queen Mary University of London.

set -ex

./fa2cmap_multi.pl -v -i tmp/contigs.fa -e BspQI -m 5 -M 80 -o tmp/

# XMAP to BED.
./bionano2Allmaps.pl -k tmp/contigs_BspQI_key.txt \
  -x tmp/BNG.xmap -o tmp/BNG_raw.bed

# List all maps that have more than four markers.
cut -f4 tmp/BNG_raw.bed | cut -d':' -f1 | sort | uniq -c \
  | awk '{if ($1 >= 4) print $2}' > tmp/BNG_informative_maps

# Grep out selected maps from BNG.bed.
grep -f tmp/BNG_informative_maps tmp/BNG_raw.bed > tmp/BNG.bed
