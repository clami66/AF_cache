#!/usr/bin/env python
import re
import os
import glob
import pickle
import string
import base64
import hashlib
import argparse
from sys import argv
from pathlib import Path
from functools import partial
from multiprocessing import Pool, TimeoutError


def convert_alignment(in_alignment, out_dir, custom_taxids=None):
    seqids = set()
    tolower = str.maketrans('', '', string.ascii_lowercase)
    print(f"Opening {in_alignment}")
    with open(in_alignment) as aln:
        a3m_data = iter(aln.readlines())
    print(f"{in_alignment} loaded")
    # always write first (target) sequence to file
    target_header = ""

    while not target_header.startswith(">"):
        target_header = next(a3m_data)
    target_seq = next(a3m_data)

    custom_taxid_dic = {}
    if custom_taxids is not None:
        with open(custom_taxids, "r") as f:
            for line in f:
                split_line = line.strip().split(" ")
                try:
                    custom_taxid_dic[int(split_line[0])] = int(split_line[1])
                except: # if the conversion is between Mnemonic IDs (e.g. 'HERPS HUMAN')
                    custom_taxid_dic[split_line[0]] = split_line[1]

    # prepare the directory structure as AF wants it
    target_id = target_header.strip().strip(">").split()[0]
    if "|" in target_id:
        target_id = target_id.split("|")[1]
    target_id = re.sub('[^0-9a-zA-Z]+', '', target_id)
    print(f"Target ID: {target_id}")
    Path(args.out_dir, target_id, "msas", "A").mkdir(parents=True, exist_ok=True)
    a3m_out = open(f"{args.out_dir}/{target_id}/msas/A/mmseqs2_hits.a3m", "w")
    a3m_out.write(f">{target_id}\n")
    a3m_out.write(target_seq)
    print(f"Starting conversion: {target_id}")
    for line in a3m_data:
        if line.startswith(">"):
            # TODO add taxid or repid conversion here
            a3m_out.write(line)
        elif line.startswith("#") or len(line.strip()) < 1:
            continue
        else:
            a3m_out.write(line)

    a3m_out.close()


def main(args):
    a3ms = glob.glob(f"{args.in_path}/*.a3m")
    print(f"Converting alignments: {a3ms}")
    with Pool(processes=int(args.n_cpu)) as pool:
        pool.map(partial(convert_alignment, out_dir=args.out_dir, custom_taxids=args.custom_taxids), a3ms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Converts an a3m alignment (e.g. MMseqs2 alignment) to a Stockholm alignment \
                                 that can be used for MSA pairing in AlphaFold-multimer runs")
    parser.add_argument("in_path", help = "Path to the a3m alignment file directory")
    parser.add_argument("out_dir", help = "Path to output directory")
    parser.add_argument("--n_cpu", default=8, type=int)
    parser.add_argument("--custom_taxids", default=None, type=str, help = "Space-separated file of taxIDs or mnemoic IDs to replace in the outputs (e.g. '12345 6789' to convert taxID 12345 to 6789), for example in case to replace virus taxIDs with their respective host taxIDs")
    args = parser.parse_args()

    main(args)
