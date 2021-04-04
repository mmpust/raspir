#!/bin/python
# Marie-Madlen Pust
# Last updated: 04.04.2021

import sys

if len(sys.argv) != 3:
    print("USAGE: emove_special_characters_in_refdb.py <input.fa> <output.fa>")
    sys.exit(1)

with open(sys.argv[1], "r") as input_fasta, open(sys.argv[2], "w") as output_fasta:
        for line in input_fasta:
                if not line.strip():
                        continue
                output_fasta.write(line.replace(".", "_").replace(" ", "_").replace(",", "_"))
