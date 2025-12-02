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
from Bio import AlignIO
import pandas as pd

seqid_to_id = dict()

def get_accession2id(db_path, taxid=True):
    if taxid:
        pkl_path = f"{db_path}/seqid2taxid.pkl"
    else:
        pkl_path = f"{db_path}/seqid2repid.pkl"

    if not os.path.exists(pkl_path):
        print("Generating accession to taxid dictionary")
        make_accession2id(db_path)
    
    print("Loading accession to taxid dictionary")
    with open(pkl_path, "rb") as input_file:
        return pickle.load(input_file)

def make_accession2id(db_path):
    print("Generating accession2taxid pickle file")
    seqid_to_taxid_dic = {}
    seqid_to_repid_dic = {}
    header_files = glob.glob(f"{db_path}/*_h.tsv")
    ptax = re.compile("TaxID=([0-9]+)")
    prep = re.compile("RepID=[A-Z0-9]+[_]([A-Z0-9]+)")

    for header_file in header_files:
        print(f"Looking for taxids in {header_file}")
        with open(header_file) as headers:
            for header in headers:
                match_taxid = ptax.search(header)
                match_repid = prep.search(header)
                accession = header.split()[1]
                if match_taxid:
                    taxid = match_taxid.group(1)
                    seqid_to_taxid_dic[accession] = int(taxid)
                if match_repid:
                    repid = match_repid.group(1)
                    seqid_to_repid_dic[accession] = repid

    with open(seqid_to_taxid, 'wb') as f:
      pickle.dump(seqid_to_taxid_dic, f, protocol=4)
    with open(seqid_to_repid, 'wb') as f:
      pickle.dump(seqid_to_repid_dic, f, protocol=4)

def convert_alignment(in_alignment, out_dir, taxid=True, duplicate=False, shortname=False, custom_taxids=None):
    seqids = set()
    tolower = str.maketrans('', '', string.ascii_lowercase)
    print(f"Opening {in_alignment}")
    with open(in_alignment) as aln:
        a3m_data = iter(aln.readlines())
    print(f"{in_alignment} loaded")
    # always write first (target) sequence to file
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
    target_id = target_header.strip().strip(">")
    print(f"Target ID: {target_id}")
    Path(args.out_dir, target_id, "msas", "A").mkdir(parents=True, exist_ok=True)
    pseudo_uniprot = open(f"{args.out_dir}/{target_id}/msas/A/uniprot_hits.a3m", "w")
    pseudo_uniprot.write(target_header)
    pseudo_uniprot.write(target_seq)
    print(f"Starting conversion: {target_id}")
    for line in a3m_data:
        if line.startswith(">"):
            memid = None
            seqid = line.split()[0].strip(">") # gets accession ID

            if seqid not in seqids:
                seqids.add(seqid)
            else:
                continue

            if seqid in seqid_to_id:
                memid = seqid_to_id[seqid]
                memid = custom_taxid_dic[memid] if memid in custom_taxid_dic else memid
                if taxid:
                    memid = base64.urlsafe_b64encode(hashlib.md5(str(memid).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:5].upper()

            if memid:
                # UniRef100_A0A2I3H6P3 -> A0A2I3H6P3
                if shortname and len(seqid) > 7:
                    seqid = re.sub("UniRef100_", "", seqid)
                    seqid = base64.urlsafe_b64encode(hashlib.md5(str(seqid).encode('utf-8')).digest()).decode("utf-8").replace("_", "").replace("-", "")[:7].lower()

                seqid = re.sub(r"[^a-zA-Z0-9]", "", seqid)
                pseudo_uniprot.write(f">tr|{seqid}|{seqid}_{memid}/1-{len(target_seq)}\n")

        elif memid:
            pseudo_uniprot.write(line.translate(tolower))

    pseudo_uniprot.close()

    input_handle  = open(f"{out_dir}/{target_id}/msas/A/uniprot_hits.a3m", "r")
    output_handle = open(f"{out_dir}/{target_id}/msas/A/uniprot_hits.sto", "w")

    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "stockholm")

    output_handle.close()
    input_handle.close()

    # bfd symlink
    src = os.path.abspath(in_alignment)
    dst = f"{os.path.abspath(out_dir)}/{target_id}/msas/A/bfd_uniref_hits.a3m"

    if os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)


def main(args):
    a3ms = glob.glob(f"{args.in_path}/*.a3m")
    print(f"Converting alignments: {a3ms}")
    with Pool(processes=int(args.n_cpu)) as pool:
        pool.map(partial(convert_alignment, out_dir=args.out_dir, taxid=taxid,
                         duplicate=args.duplicate, shortname=args.shortname, custom_taxids=args.custom_taxids), a3ms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Converts an a3m alignment (e.g. MMseqs2 alignment) to a Stockholm alignment \
                                 that can be used for MSA pairing in AlphaFold-multimer runs")
    parser.add_argument("in_path", help = "Path to the a3m alignment file directory")
    parser.add_argument("out_dir", help = "Path to output directory")
    parser.add_argument("db_path", help = "Path to ColabFold databases")
    parser.add_argument("--n_cpu", default=8, type=int)
    parser.add_argument("--taxid", action="store_true", default=True, help = "If tax IDs should be used in the pairing procedure")
    parser.add_argument("--repid", action="store_true", default=False, help = "If mnemonic species IDs should be used in the pairing procedure")
    parser.add_argument("--duplicate", action="store_true", default=False, help = "If we should add duplicate species to the uniprot file")
    parser.add_argument("--shortname", action="store_true", default=False, help = "To replace long names in the database with 7-char names - necessary for pairing")
    parser.add_argument("--custom_taxids", default=None, type=str, help = "Space-separated file of taxIDs or mnemoic IDs to replace in the outputs (e.g. '12345 6789' to convert taxID 12345 to 6789), for example in case to replace virus taxIDs with their respective host taxIDs")
    args = parser.parse_args()
    taxid = False if args.repid else True
    seqid_to_id = get_accession2id(args.db_path, taxid=taxid)
    main(args)
