#!/usr/bin/env python

import os
import sys
import math
import glob
import argparse
import subprocess
from pathlib import Path
from shutil import which, move, rmtree
from alphafold.data.tools import hhblits, jackhmmer, mmseqs2
from alphafold.data.pipeline import run_msa_tool
import numpy as np

parser = argparse.ArgumentParser(
    description="Run one of the MSA tools used in AlphaFold. This includes mmseqs2 alignments as run in ColabFold"
)
parser.add_argument("fasta_path", help="Path to the fasta file for a single protein")
parser.add_argument(
    "alignment_type",
    choices=[
        "uniref90",
        "mgnify",
        "uniref30",
        "bfd",
        "bfd_small",
        "uniprot",
        "mmseqs2",
    ],
)
parser.add_argument("db_path", help="Path to database for alignment")
parser.add_argument("--out_dir", type=str, help="Output path for the MSA file")
parser.add_argument("--max_hits", default=None, help="Number of max hits to report")
parser.add_argument("--use_precomputed_msas", action="store_true")
parser.add_argument("--n_cpu", default=32, type=int)
parser.add_argument("--n_iter", type=int, default=None)
parser.add_argument("--mmseqs", default="mmseqs", help="mmseqs binary path")
parser.add_argument(
    "--db1", default="uniref30_2302_db", help="ColabFold uniref DB name for MMseqs2"
)
parser.add_argument(
    "--db2",
    default="colabfold_envdb_202108_db",
    help="ColabFold env DB name for MMseqs2",
)
parser.add_argument("--use-env", action="store_true", help="Use env DB in MMseqs2")
parser.add_argument("--gpu", action="store_true", help="Use GPU-accelerated MMseqs2")
parser.add_argument(
    "--gpu_server", action="store_true", help="Launch MMseqs2 GPU server"
)
parser.add_argument(
    "--max_accept", default=100000, help="MMseqs2 max accepted alignments"
)

args = parser.parse_args()

fasta_path = args.fasta_path
db_path = args.db_path
alignment_type = args.alignment_type
msa_out_dir = args.out_dir
n_cpu = args.n_cpu
n_iter = args.n_iter

if not args.max_hits:
    if alignment_type == "mgnify":
        max_hits = 501
    elif alignment_type in ["uniprot", "mmseqs", "mmseqs2"]:
        max_hits = None
    else:
        max_hits = 10000
else:
    max_hits = args.max_hits

if alignment_type in ["bfd_small", "uniref90", "mgnify", "uniprot"]:
    if not n_iter:
        n_iter = 1
    runner = jackhmmer.Jackhmmer(
        binary_path=which("jackhmmer"),
        database_path=db_path,
        n_cpu=n_cpu,
        n_iter=n_iter,
    )
elif alignment_type in ["uniref30", "bfd"]:
    if not n_iter:
        n_iter = 3
    runner = hhblits.HHBlits(
        binary_path=which("hhblits"), databases=[db_path], n_cpu=n_cpu, n_iter=n_iter
    )
elif alignment_type in ["mmseqs", "mmseqs2"]:
    Path(msa_out_dir).mkdir(parents=True, exist_ok=True)
    runner = mmseqs2.MMseqs2(
        binary_path=Path(args.mmseqs) if args.mmseqs else which("mmseqs"),
        uniref_db=db_path + args.db1,
        metagenomic_db=db_path + args.db2 if args.use_env else None,
        msa_out_dir=Path(msa_out_dir),
        n_cpu=n_cpu,
        gpu=args.gpu,
        gpu_server=args.gpu_server,
        max_accept=args.max_accept,
    )


msa_format = (
    "a3m" if alignment_type in ["uniref30", "bfd", "mmseqs", "mmseqs2"] else "sto"
)
msa_out_file = f"{alignment_type}_hits.{msa_format}"
msa_out_path = os.path.join(msa_out_dir, msa_out_file)

_ = run_msa_tool(
    msa_runner=runner,
    input_fasta_path=fasta_path,
    msa_out_path=msa_out_path,
    msa_format=msa_format,
    use_precomputed_msas=args.use_precomputed_msas,
    max_sto_sequences=max_hits,
)

if alignment_type in ["mmseqs", "mmseqs2"]:
    os.remove(msa_out_path)
    if os.path.exists(f"{msa_out_dir}/alignments"):
        rmtree(f"{msa_out_dir}/alignments")
    alignments_dir = glob.glob(f"{msa_out_dir}/alignments?*")[
        0
    ]  # only one temp output dir should be there
    move(
        alignments_dir, f"{msa_out_dir}/alignments"
    )  # older alignments will be overwritten
