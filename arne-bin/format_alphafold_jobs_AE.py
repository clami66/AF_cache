import os
import argparse
from sys import argv
from math import inf
from glob import glob
from pathlib import Path
from bisect import bisect
from itertools import combinations_with_replacement,combinations,product

from Bio import SeqIO

def get_slurm_profile(proj_id, max_len, log_path):

    return f"""#!/bin/bash
#SBATCH -A {proj_id}
#SBATCH --gpus 1
#SBATCH -t 3-00:00:00
#SBATCH -C {'thin' if max_len != inf else 'fat'}
#SBATCH -o {log_path}

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='6.0'

module load Anaconda/2021.05-nsc1
"""

def get_file_list(file):
    f = open(file, 'r')
    pairlist=[]
    tmplist=f.readlines()
    for i in tmplist:
        pairlist.append(i.strip("\n"))
    f.close()
    return pairlist

def get_fasta_records(fasta_files: list):
    fasta_records = []
    for ff in fasta_files:
        fasta_path = Path(ff)
        with open(fasta_path) as f:
            record = next(SeqIO.parse(f, "fasta"))
            fasta_records.append((record, fasta_path.stem))
    return fasta_records


def estimate_gpu_runtime(seqlen, lam=0.0001): # this is a rough estimate and largely dependent on the GPU speed/VRAM size
    return seqlen**2 * lam / 3600


def format_af_command(target_list, out_dir, pad_to_size=None, pickle_dir=None, flagfile=f"/proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/multimer_all_vs_all.flag",other_args=""):
    scripts_path = os.path.dirname(os.path.realpath(__file__))
    #flagfile = f"/proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/multimer_all_vs_all.flag"
    pickle_flag = f"--pickle_cache {pickle_dir}" if pickle_dir else ""
    pad_flag = f"--pad_to_size {pad_to_size}" if pad_to_size else ""
    condaenv = "AF_cache" if pad_to_size else "af_server"
    #return f"conda activate /proj/beyondfold/apps/.conda/envs/af_server \n python /proj/beyondfold/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile {flagfile} --output_dir {out_dir} --fasta_paths {','.join(target_list)} {pickle_flag} {pad_flag} {' '.join(other_args)}"
    return f"conda activate /proj/beyondfold/apps/.conda/envs/af_server \n python /proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile {flagfile} --output_dir {out_dir} --fasta_paths {','.join(target_list)} {pickle_flag} {pad_flag} {' '.join(other_args)}"



def define_pairs(fasta_records, out_dir, splits, list, write_fastas=False, overwrite_output=True,include_homomers=True,both_directions=False):


    if (len(list)>0 or both_directions):
        all_pairs = product(fasta_records, repeat=2)
    elif(include_homomers):
        all_pairs = combinations_with_replacement(fasta_records, 2)
    else:
        all_pairs = combinations(fasta_records, 2)

    pair_bins = {split:[] for split in splits}
    for pair in all_pairs:
        pair_id = f"{pair[0][1]}_{pair[1][1]}"
        pair_records = (pair[0][0], pair[1][0])
        if (len(list)==0 or pair_id in list) :
            af_output = glob(f"{out_dir}/{pair_id}/unrelaxed*pdb")

            if not af_output or overwrite_output:
                pair_fasta = Path(out_dir, pair_id, f"{pair_id}.fasta")
                if args.write_fastas:
                    pair_folder = Path(out_dir, pair_id)
                    pair_folder.mkdir(parents=True,exist_ok=True)
                        
                    with open(pair_fasta, "w") as pf:
                        SeqIO.write(pair_records[0], pf, "fasta")
                        SeqIO.write(pair_records[1], pf, "fasta")

                pair_size = len(pair_records[0].seq) + len(pair_records[1].seq)

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

    fasta_records = get_fasta_records(glob(f"{args.in_path}/*.fasta"))
    #if (args.file_list):
    try:
        pairlist=get_file_list(args.file_list)
    except:
        pairlist=[]
    binned_pairs = define_pairs(fasta_records, out_dir, splits,pairlist,write_fastas=args.write_fastas, overwrite_output=args.overwrite_output,
                              include_homomers=args.include_homomers,both_directions=args.both_directions)
    
    Path(out_dir, "sbatch_scripts").mkdir(parents=True, exist_ok=True)
    Path(out_dir, "logs").mkdir(parents=True, exist_ok=True)
    estimated_gpu_runtime = 0
    num_jobs = 0
    for (max_len, this_bin), max_size in zip(binned_pairs.items(), max_job_size):
        num_targets = len(this_bin)
        for chunk_n, index in enumerate(range(0, num_targets,max_size )):
            target_chunk = this_bin[index:index + max_size]
            target_fastas = [target[0] for target in target_chunk]
            target_sizes = [target[1] for target in target_chunk]

            pad_to_size = f"{max_len},{max_depth}" if max_size > 1 else f"{target_sizes[0]},{max_depth}"

            num_jobs += 1
            estimated_gpu_runtime += 120/3600 # ~2 minutes per job to compile models etc
            if not args.estimate_gpu_runtime:
                command_file = Path(out_dir, "sbatch_scripts", f"{max_len}_{chunk_n}.sh")
                log_file = Path(out_dir, "logs", f"{max_len}_{chunk_n}.log")
                with open(command_file, "w") as command:
                    command.write(get_slurm_profile(args.proj_id, max_len, str(log_file)))
                    command.write("\n")
                    command.write(format_af_command([target[0] for target in target_chunk], out_dir, pickle_dir=args.pickle_dir, pad_to_size=pad_to_size, flagfile=args.flagfile,other_args=af_args))
                    command.write("\n")

    if args.estimate_gpu_runtime:
        print(f"(Under)estimated GPU core hours: {int(estimated_gpu_runtime)} core hours for {num_jobs} jobs")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Format all vs. all AlphaFold job commands given a set of fasta files")
    parser.add_argument("in_path", help = "Path to the directory containing the fasta files")
    parser.add_argument("--file_list", help = "Path to file containing a list of files to run (if not desire all against all)",default="")
    parser.add_argument("--include_homomers", action="store_true", default=False, help="Also include homomers")
    parser.add_argument("--both_directions", action="store_true", default=False, help="Run AB as well as BA")
    parser.add_argument("out_dir", help = "Path to output directory (as will be used in AlphaFold)")
    parser.add_argument("--flagfile", help = "Flagfile with parameters to AF", default=f"/proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/multimer_all_vs_all.flag")
    parser.add_argument("--pickle_dir", default="", help="Path to directory containing pickled features for all monomers in set")
    parser.add_argument("--proj_id", default="berzelius-2023-328", help="SLURM project ID")
    parser.add_argument("--write_fastas", action="store_true", default=False, help="If the fasta files and folder structure for all pairs should be initialized")
    parser.add_argument("--overwrite_output", action="store_true", default=False, help="If previously generated dimer predictions should be overwritten")
    parser.add_argument("--splits", nargs="+", default=[400, 800, 1000, 1200, 1400, 1600, 4500], help="Boundaries (sum of sequences length) to group multiple inference jobs")
    parser.add_argument("--max_job_size", nargs="+", default=[1000, 500, 100, 100, 100, 50, 1], help="When grouping jobs by length (with --splits), max number of targets that should run on the same AF python command for each split")
    parser.add_argument("--estimate_gpu_runtime", action="store_true")

    args, unknownargs = parser.parse_known_args()

    main(args, unknownargs)
