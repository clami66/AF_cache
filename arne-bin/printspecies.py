#!/usr/bin/env python3

import re
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import string
import base64
import hashlib
from sys import argv
import glob
from pathlib import Path


parser = argparse.ArgumentParser(description="Replace Host ID with Virus IDS")
parser.add_argument("--infile", help = "Path to the taxidfile",default="/proj/beyondfold/apps/colabfold_databases/seqid2taxid.pkl")
parser.add_argument("--msa",help="MSAdir",required=True)
parser.add_argument("--out",help="Output file",required=True)
args, unknownargs = parser.parse_known_args()

with open(args.infile, "rb") as file:
    seqid_to_id = pickle.load(file)


a3ms = glob.glob(f"{args.msa}/*.a3m")
seqids = set()
tolower = str.maketrans('', '', string.ascii_lowercase)
for in_alignment in a3ms:
    print(f"Opening {in_alignment}")
    with open(in_alignment) as aln:
        a3m_data = iter(aln.readlines())
    print(f"{in_alignment} loaded")
    target_header = next(a3m_data)
    target_seq = next(a3m_data)
    
    target_id = target_header.strip().strip(">")
    print(f"Target ID: {target_id}")
    Path(args.out).mkdir(parents=True, exist_ok=True)
    pseudo_uniprot = open(f"{args.out}/{target_id}.a3m", "w")
    pseudo_uniprot.write(target_header)
    pseudo_uniprot.write(target_seq)
    print(f"Starting conversion: {target_id}")

    for line in a3m_data:
        if line.startswith(">"):
            memid = None
            seqid = line.split()[0].strip(">") # gets accession ID

            if seqid not in seqids: # avoids duplicates
                seqids.add(seqid)
            else:
                continue
            
            if seqid in seqid_to_id:
                memid = seqid_to_id[seqid]

            taxid = base64.urlsafe_b64encode(hashlib.md5(str(memid).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:5].upper()
                
            if memid:
                seqid_alpha = re.sub(r"[^a-zA-Z0-9]", '', seqid)
                pseudo_uniprot.write(f">tr|{seqid_alpha}|{seqid_alpha}_{taxid}/1-{len(target_seq)} {memid}\n")
        elif memid:
            pseudo_uniprot.write(line.translate(tolower))
