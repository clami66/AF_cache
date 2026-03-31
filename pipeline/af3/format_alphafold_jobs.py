import os
import json
import hashlib
import argparse
from sys import argv, exit
from math import inf
from glob import glob
from pathlib import Path
from bisect import bisect
from itertools import combinations_with_replacement, combinations, product

from Bio import SeqIO
from string import ascii_uppercase, ascii_lowercase
ascii_upperlower = ascii_uppercase + ascii_lowercase

def bash_header():

    return f"""#!/bin/bash
module load Mambaforge/23.3.1-1-hpc1-bdist

"""


def get_fasta_record(fasta_file):
    try:
        with open(fasta_file) as f:
            record = next(SeqIO.parse(f, "fasta"))
            sequence_hash = hashlib.md5(record.seq.encode()).hexdigest()
    except:
        print(f"Missing fasta record: {fasta_file}")
        return None
    return record, sequence_hash


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


def estimate_gpu_runtime(seqlen, lam=0.0001): # this is a rough estimate and largely dependent on the GPU speed/VRAM size
    return seqlen**2 * lam / 3600


def format_af_command(json_input_dir, out_dir, flagfiles=None, af3_path="alphafold3/", other_args=""):
    flag_param = f"--flagfile {' --flagfile '.join(flagfiles)}" if flagfiles is not None else ""
    return f"python {af3_path}/run_alphafold.py --output_dir {out_dir} --input_dir {json_input_dir} {flag_param} {' '.join(other_args)}"


def merge_jsons(jsons):
    parsed_jsons = [json.load(open(j)) for j in jsons]

    merged_json = parsed_jsons[0]

    for i, extra_json in enumerate(parsed_jsons[1:]):
        chain = ascii_upperlower[i+1]
        sequences = extra_json["sequences"]
        sequences[0]["protein"]["id"] = chain
        sequences[0]["protein"]["pairedMsa"] = ""
        merged_json["sequences"].append(sequences[0])

    return merged_json


def group_multimers(fasta_records, out_dir, splits, multimer_list, json_dir, write_fastas=False, overwrite_output=True, include_homomers=True, both_directions=False):

    if multimer_list:
        all_multimers = multimer_list
    elif both_directions:
        all_multimers = product(fasta_records, repeat=2)
    elif include_homomers:
        all_multimers = combinations_with_replacement(fasta_records, 2)
    else:
        all_multimers = combinations(fasta_records, 2)

    multimer_bins = {split:[] for split in splits}
    multimer_json_dir = Path(out_dir, "multimer_jsons")
    multimer_json_dir.mkdir(parents=True, exist_ok=True)
    for count, multimer in enumerate(all_multimers):
        if not count % 1000:
            print(count)
        multimer_id = "_".join([p[1] for p in multimer])
        multimer_records = [p[0][0] for p in multimer]
        multimer_hashes = [p[0][1] for p in multimer]
        merged_json = merge_jsons([f"{json_dir}/{h}.json" for h in multimer_hashes])
        merged_json["name"] = multimer_id


        multimer_json = Path(multimer_json_dir, f"{multimer_id}.json")
        with open(multimer_json, "w") as out:
            json.dump(merged_json, out, indent=4, sort_keys=True)

        af_output = glob(f"{out_dir}/{multimer_id}/unrelaxed*pdb")

        if not af_output or overwrite_output:
            multimer_fasta = Path(out_dir, multimer_id, f"{multimer_id}.fasta")
            if args.write_fastas:
                with open(multimer_fasta, "w") as pf:
                    for pr in multimer_records:
                        SeqIO.write(pr, pf, "fasta")

            multimer_size = sum([len(pr.seq) for pr in multimer_records])

            multimer_bin = splits[bisect(splits, multimer_size)]
            multimer_bins[multimer_bin].append((str(multimer_json), multimer_size))
    return multimer_bins


def main(args, af_args):
    out_dir = str(Path(args.out_dir).resolve())
    splits = [int(split) for split in args.splits]
    splits.append(inf)
    max_job_size = [int(jobs) for jobs in args.max_job_size]
    max_job_size.append(1)

    if args.file_list:
        multimer_list = get_records_from_list(args.file_list, args.in_path, sep=args.list_separator)
        fasta_records = []
    else:
        fasta_records = get_records_from_dir(glob(f"{args.in_path}/*.fasta"))
        multimer_list = []

    multimers = group_multimers(fasta_records, out_dir, splits, multimer_list, json_dir=args.json_dir, write_fastas=args.write_fastas, overwrite_output=args.overwrite_output, include_homomers=args.include_homomers, both_directions=args.both_directions)
    Path(out_dir, "sbatch_scripts").mkdir(parents=True, exist_ok=True)
    Path(out_dir, "logs").mkdir(parents=True, exist_ok=True)

    num_jobs = 0
    for (max_len, this_bin), max_size in zip(multimers.items(), max_job_size):
        num_targets = len(this_bin)
        for chunk_n, index in enumerate(range(0, num_targets, max_size)):
            target_chunk = this_bin[index:index + max_size]
            target_jsons = [target[0] for target in target_chunk]
            target_sizes = [target[1] for target in target_chunk]

            # a directory of jsons to be passed to an AF3 call so that the proteins are run in the same experiment
            input_json_dir = Path(out_dir, "bucketed_jsons", f"{max_len}_{chunk_n}")
            input_json_dir.mkdir(parents=True, exist_ok=True)

            for target_json in target_jsons:
                os.symlink(target_json, f"{input_json_dir}/{os.path.basename(target_json)}")
            num_jobs += 1

            command_file = Path(out_dir, "sbatch_scripts", f"{max_len}_{chunk_n}.sh")
            log_file = Path(out_dir, "logs", f"{max_len}_{chunk_n}.log")
            with open(command_file, "w") as command:
                command.write(bash_header())
                command.write(format_af_command(str(input_json_dir), f"{out_dir}/{max_len}_{chunk_n}", flagfiles=args.flagfiles, af3_path=args.af3_path, other_args=af_args))
                command.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Format all vs. all AlphaFold job commands given a set of fasta files")
    parser.add_argument("in_path", help = "Path to the directory containing the fasta files")
    parser.add_argument("out_dir", help = "Path to output directory (as will be used in AlphaFold)")
    parser.add_argument("--af3_path", default="alphafold3/", help="AlphaFold install path")
    parser.add_argument("--flagfiles", nargs="+", help = "Flagfile with parameters to AF", default=None)
    parser.add_argument("--file_list", help = "Path to file containing a list of files to run (if not desire all against all)",default="")
    parser.add_argument("--list_separator", default=" ", help="Character used to separate protein pairs in file list (--file_list)")
    parser.add_argument("--include_homomers", action="store_true", default=False, help="Also include homomers")
    parser.add_argument("--both_directions", action="store_true", default=False, help="Run AB as well as BA")
    parser.add_argument("--json_dir", default="", help="Path to directory containing json features for all monomers in set")
    parser.add_argument("--write_fastas", action="store_true", default=False, help="If the fasta files and folder structure for all pairs should be initialized")
    parser.add_argument("--overwrite_output", action="store_true", default=False, help="If previously generated dimer predictions should be overwritten")
    parser.add_argument("--splits", nargs="+", default=[256, 512, 768, 1024, 1280, 1536, 2048, 2560, 3072, 3584, 4096, 4608, 5120], help="Bucket boundaries to group multiple inference jobs")
    parser.add_argument("--max_job_size", nargs="+", default=[1000, 500, 100, 100, 100, 50, 1, 1, 1, 1, 1, 1, 1], help="When grouping jobs by length (with --splits), max number of targets that should run on the same AF python command for each split")

    args, unknownargs = parser.parse_known_args()

    main(args, unknownargs)
