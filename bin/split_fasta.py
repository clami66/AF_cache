#!/usr/bin/env python

import re
from sys import argv

fasta_path = argv[1]
out_path = argv[2]
target_id = ""
fid = None

with open(fasta_path) as fasta_in:
    for line in fasta_in:
        if line.startswith(">"):  # next sequence, write rest of fasta
            if fid is not None:
                fid.close()
            target_id = line.split()[0][1:]
            if "|" in target_id:
                target_id = target_id.split("|")[1]
            target_id = re.sub("[^0-9a-zA-Z]+", "", target_id)
            fid = open(f"{out_path}/{target_id}.fasta", "w")
            fid.write(f"{line}")

        else:
            fid.write(line)
