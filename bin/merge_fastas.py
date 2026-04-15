#!/usr/bin/env python
import glob
import argparse
from pathlib import Path
from Bio import SeqIO


def get_fasta_record(fasta_file):
    try:
        with open(fasta_file) as f:
            record = next(SeqIO.parse(f, "fasta"))
    except:
        print(f"Missing fasta record: {fasta_file}")
        return None
    return record


def main(args):
    fasta_list = glob.glob(f"{args.fasta}/*.fasta")

    fasta_records = [get_fasta_record(f).seq for f in fasta_list]
    fasta_ids = [Path(f).stem for f in fasta_list]

    seqs = []
    ids = []
    for fr, id in zip(fasta_records, fasta_ids):
        if str(fr) not in seqs or not args.unique:
            seqs.append(str(fr))
            ids.append(id)
        else:
            print(f"Warning: skipping duplicate sequence {str(fr)} with id {id}")

    with open(args.out_fasta, "w") as pf:
        for fr, id in zip(seqs, ids):
            pf.write(f">{id}\n")
            pf.write(f"{fr}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merges multiple inputs to fasta files into a single"
    )
    parser.add_argument("fasta", help="Directory containing multiple fasta files")
    parser.add_argument("out_fasta", help="Path to merged fasta file")
    parser.add_argument(
        "--unique",
        action="store_true",
        help="Remove duplicate sequences in merged file",
    )
    args = parser.parse_args()

    main(args)
