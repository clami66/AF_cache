#!/usr/bin/env python

import os
import argparse
from sys import argv, exit
from math import inf
from glob import glob
from pathlib import Path
from bisect import bisect
from itertools import combinations_with_replacement, combinations, product

from Bio import SeqIO


def get_fasta_record(fasta_file):
    try:
        with open(fasta_file) as f:
            record = next(SeqIO.parse(f, "fasta"))
    except:
        print(f"Missing fasta record: {fasta_file}")
        return None
    return record


def get_records_from_list(list_file, fasta_path, sep):
    with open(list_file, "r") as f:
        tmplist = f.readlines()
    multimer_list = []
    for i in tmplist:
        monomers = i.strip("\n").split(sep)
        fastas = [Path(fasta_path, f"{p}.fasta") for p in monomers]
        records = [get_fasta_record(fasta) for fasta in fastas]
        multimer_list.append([(r, p) for r, p in zip(records, monomers)])
    return multimer_list


def get_records_from_dir(fasta_files: list):
    fasta_records = []
    for ff in fasta_files:
        fasta_path = Path(ff)
        record = get_fasta_record(fasta_path)

        if record:
            fasta_records.append((record, fasta_path.stem))
    return fasta_records


def define_pairs(
    fasta_records,
    out_dir,
    splits,
    pair_list,
    write_fastas=False,
    overwrite_output=True,
    include_homomers=True,
    both_directions=False,
):
    if pair_list:
        all_pairs = pair_list
    elif both_directions:
        all_pairs = product(fasta_records, repeat=2)
    elif include_homomers:
        all_pairs = combinations_with_replacement(fasta_records, 2)
    else:
        all_pairs = combinations(fasta_records, 2)

    pair_bins = {split: [] for split in splits}
    for count, pair in enumerate(all_pairs):
        if not count % 1000:
            print(count)

        pair_id = "_".join([p[1] for p in pair])
        pair_records = [p[0] for p in pair]
        af_output = glob(f"{out_dir}/{pair_id}/unrelaxed*pdb")

        if not af_output or overwrite_output:
            pair_fasta = Path(out_dir, pair_id, f"{pair_id}.fasta")
            if write_fastas:
                pair_folder = Path(out_dir, pair_id)
                pair_folder.mkdir(parents=True, exist_ok=True)

                with open(pair_fasta, "w") as pf:
                    for pr in pair_records:
                        SeqIO.write(pr, pf, "fasta")

            pair_size = sum([len(pr.seq) for pr in pair_records])

            pair_bin = splits[bisect(splits, pair_size)]
            pair_bins[pair_bin].append((str(pair_fasta), pair_size))
    return pair_bins


def main(args, af_args):
    splits = [int(split) for split in args.splits]
    splits.append(inf)
    max_job_size = [int(jobs) for jobs in args.max_job_size]
    max_job_size.append(1)
    max_depth = 20000

    out_dir = str(Path(args.out_dir).resolve())

    if args.file_list:
        pairlist = get_records_from_list(
            args.file_list, args.in_path, sep=args.list_separator
        )
        fasta_records = []
    else:
        fasta_records = get_records_from_dir(glob(f"{args.in_path}/*.fasta"))
        pairlist = []

    binned_pairs = define_pairs(
        fasta_records,
        out_dir,
        splits,
        pairlist,
        write_fastas=args.write_fastas,
        overwrite_output=args.overwrite_output,
        include_homomers=args.include_homomers,
        both_directions=args.both_directions,
    )

    Path(out_dir, "logs").mkdir(parents=True, exist_ok=True)
    for (max_len, this_bin), max_size in zip(binned_pairs.items(), max_job_size):
        num_targets = len(this_bin)
        for chunk_n, index in enumerate(range(0, num_targets, max_size)):
            target_chunk = this_bin[index : index + max_size]
            target_fastas = [target[0] for target in target_chunk]
            target_sizes = [target[1] for target in target_chunk]

            pad_to_size = (
                f"{max_len},{max_depth}"
                if max_size > 1
                else f"{target_sizes[0]},{max_depth}"
            )

            input_fasta_dir = Path(f"chunk_{max_len}_{chunk_n}")
            input_fasta_dir.mkdir(parents=True, exist_ok=True)

            flag_file = Path(input_fasta_dir, "chunk.flags")
            other_flags = Path(input_fasta_dir, "other.flags")

            with open(flag_file, "w") as flags:
                fasta_inputs = [
                    f"{input_fasta_dir}/{os.path.basename(tf)}" for tf in target_fastas
                ]
                flags.write(f"--fasta_paths={','.join(fasta_inputs)}\n")
                for target_fasta, fasta_input in zip(target_fastas, fasta_inputs):
                    os.symlink(target_fasta, fasta_input)
                flags.write(f"--pad_to_size={pad_to_size}\n")
                flags.write(f"--pickle_cache={args.pickle_dir}\n")
                flags.write(f"--flagfile={args.flagfile}\n")

            with open(other_flags, "w") as flags:
                flags.write(f"{' '.join(af_args)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Format all vs. all AlphaFold job commands given a set of fasta files"
    )
    parser.add_argument(
        "in_path", help="Path to the directory containing the fasta files"
    )
    parser.add_argument(
        "out_dir", help="Path to output directory (as will be used in AlphaFold)"
    )
    parser.add_argument(
        "--file_list",
        help="Path to file containing a list of files to run (if not desire all against all)",
        default="",
    )
    parser.add_argument(
        "--list_separator",
        default=" ",
        help="Character used to separate protein pairs in file list (--file_list)",
    )
    parser.add_argument(
        "--include_homomers",
        action="store_true",
        default=False,
        help="Also include homomers",
    )
    parser.add_argument(
        "--both_directions",
        action="store_true",
        default=False,
        help="Run AB as well as BA",
    )
    parser.add_argument(
        "--flagfile", help="Flagfile with parameters to AF", default=None
    )
    parser.add_argument(
        "--pickle_dir",
        default="",
        help="Path to directory containing pickled features for all monomers in set",
    )
    parser.add_argument(
        "--write_fastas",
        action="store_true",
        default=False,
        help="If the fasta files and folder structure for all pairs should be initialized",
    )
    parser.add_argument(
        "--overwrite_output",
        action="store_true",
        default=False,
        help="If previously generated dimer predictions should be overwritten",
    )
    parser.add_argument(
        "--splits",
        nargs="+",
        default=[400, 800, 1000, 1200, 1400, 1600, 4500],
        help="Boundaries (sum of sequences length) to group multiple inference jobs",
    )
    parser.add_argument(
        "--max_job_size",
        nargs="+",
        default=[1000, 500, 100, 100, 100, 50, 1],
        help="When grouping jobs by length (with --splits), max number of targets that should run on the same AF python command for each split",
    )

    args, unknownargs = parser.parse_known_args()

    main(args, unknownargs)
