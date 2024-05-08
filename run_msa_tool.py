import os
import sys
import math
import argparse
import subprocess
from pathlib import Path
from shutil import which, rmtree
from absl import logging
from typing import List, Union
from alphafold.data.tools import hhblits
from alphafold.data.tools import jackhmmer
from alphafold.data.pipeline import run_msa_tool
import numpy as np

logging.set_verbosity(logging.INFO)
logging.use_absl_handler()

def run_mmseqs(mmseqs: Path, params: List[Union[str, Path]]):
    params_log = " ".join(str(i) for i in params)
    logging.info(f"Running {mmseqs} {params_log}")
    subprocess.check_call([mmseqs] + params)

def mmseqs_search_monomer(
    dbbase: Path,
    base: Path,
    uniref_db: Path = Path("uniref30_2202_db"),
    metagenomic_db: Path = Path("colabfold_envdb_202108_db"),
    mmseqs: Path = Path("mmseqs"),
    use_env: bool = True,
    use_templates: bool = False,
    filter: bool = True,
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    s: float = 8,
    db_load_mode: int = 2,
    threads: int = 32,
    split_memory_limit: str = "",
):
    """Run mmseqs with a local colabfold database set

    db1: uniprot db (UniRef30)
    db2: Template (unused by default)
    db3: metagenomic db (colabfold_envdb_202108 or bfd_mgy_colabfold, the former is preferred)
    """
    if filter:
        # 0.1 was not used in benchmarks due to POSIX shell bug in line above
        #  EXPAND_EVAL=0.1
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    used_dbs = [uniref_db]
    if use_env:
        used_dbs.append(metagenomic_db)

    for db in used_dbs:
        if not dbbase.joinpath(f"{db}.dbtype").is_file():
            raise FileNotFoundError(f"Database {db} does not exist")
        if (
            not dbbase.joinpath(f"{db}.idx").is_file()
            and not dbbase.joinpath(f"{db}.idx.index").is_file()
        ):
            logger.info("Search does not use index")
            db_load_mode = 0
            dbSuffix1 = "_seq"
            dbSuffix2 = "_aln"
        else:
            dbSuffix1 = ".idx"
            dbSuffix2 = ".idx"

    # fmt: off
    # @formatter:off
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-s", str(s), "-e", "0.1", "--max-seqs", "10000",]
    if split_memory_limit:
        search_param.extend(["--split-memory-limit", split_memory_limit, "--split", "0"])
    filter_param = ["--filter-msa", str(filter), "--filter-min-enable", "1000", "--diff", str(diff), "--qid", "0.0,0.2,0.4,0.6,0.8,1.0", "--qsc", "0", "--max-seq-id", "0.95",]
    expand_param = ["--expansion-mode", "0", "-e", str(expand_eval), "--expand-filter-clusters", str(filter), "--max-seq-id", "0.95",]

    run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), dbbase.joinpath(uniref_db), base.joinpath("res"), base.joinpath("tmp"), "--threads", str(threads)] + search_param)
    run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res"), dbbase.joinpath(f"{uniref_db}{dbSuffix2}"), base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + expand_param)
    run_mmseqs(mmseqs, ["mvdb", base.joinpath("tmp/latest/profile_1"), base.joinpath("prof_res")])
    run_mmseqs(mmseqs, ["lndb", base.joinpath("qdb_h"), base.joinpath("prof_res_h")])
    run_mmseqs(mmseqs, ["align", base.joinpath("prof_res"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp"), base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads", str(threads), "--alt-ali", "10", "-a"])
    run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"),
                        base.joinpath("res_exp_realign"), base.joinpath("res_exp_realign_filter"), "--db-load-mode",
                        str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0", "--threads",
                        str(threads), "--max-seq-id", "1.0", "--filter-min-enable", "100"])
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"),
                        base.joinpath("res_exp_realign_filter"), base.joinpath("uniref.a3m"), "--msa-format-mode",
                        "6", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign_filter")])

    if use_env:
        run_mmseqs(mmseqs, ["search", base.joinpath("prof_res"), dbbase.joinpath(metagenomic_db), base.joinpath("res_env"), base.joinpath("tmp"), "--threads", str(threads)] + search_param)
        run_mmseqs(mmseqs, ["expandaln", base.joinpath("prof_res"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"), base.joinpath("res_env"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix2}"), base.joinpath("res_env_exp"), "-e", str(expand_eval), "--expansion-mode", "0", "--db-load-mode", str(db_load_mode), "--threads", str(threads)])
        run_mmseqs(mmseqs,
                   ["align", base.joinpath("tmp/latest/profile_1"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                    base.joinpath("res_env_exp"), base.joinpath("res_env_exp_realign"), "--db-load-mode",
                    str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads",
                    str(threads), "--alt-ali", "10", "-a"])
        run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp_realign"), base.joinpath("res_env_exp_realign_filter"),
                            "--db-load-mode", str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0",
                            "--max-seq-id", "1.0", "--threads", str(threads), "--filter-min-enable", "100"])
        run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp_realign_filter"),
                            base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m"), "--msa-format-mode", "6",
                            "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)


        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env_exp_realign_filter")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env_exp_realign")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env_exp")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env")])

    if use_env:
        run_mmseqs(mmseqs, ["mergedbs", base.joinpath("qdb"), base.joinpath("final.a3m"), base.joinpath("uniref.a3m"), base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
    else:
        run_mmseqs(mmseqs, ["mvdb", base.joinpath("uniref.a3m"), base.joinpath("final.a3m")])

    run_mmseqs(mmseqs, ["unpackdb", base.joinpath("final.a3m"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".a3m"])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("final.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("uniref.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res")])

    for file in base.glob("prof_res*"):
        file.unlink()
    rmtree(base.joinpath("tmp"))

parser = argparse.ArgumentParser(description="Run one of the MSA tools used in AlphaFold. This includes mmseqs2 alignments as run in ColabFold")
parser.add_argument("fasta_path", help = "Path to the fasta file for a single protein")
parser.add_argument("alignment_type", choices=["uniref90", "mgnify", "uniref30", "bfd", "bfd_small", "uniprot", "mmseqs2"])
parser.add_argument("db_path", help = "Path to database for alignment")
parser.add_argument("--out_dir", type = str, help = "Output path for the MSA file")
parser.add_argument("--max_hits", default=None, help="Number of max hits to report")
parser.add_argument("--use_precomputed_msas", action="store_true")
parser.add_argument("--n_cpu", default=8)
parser.add_argument("--n_iter", type = int, default=None)
parser.add_argument("--colab_dir", default="/proj/beyondfold/apps/ColabFold/", help="ColabFold code directory to run MMseqs2 alignments")
parser.add_argument("--mmseqs", default="mmseqs", help="mmseqs binary path")
parser.add_argument("--db1", type=Path, default=Path("uniref30_2302_db"), help="ColabFold uniref DB name for MMseqs2")
parser.add_argument("--db2", type=Path, default=Path("colabfold_envdb_202108_db"), help="ColabFold env DB name for MMseqs2")
parser.add_argument("--use-env", type=int, default=1, choices=[0, 1], help="Use env DB in MMseqs2")
parser.add_argument("--db-load-mode", type=int, default=0, help="DB loading mode in MMseqs2")
parser.add_argument("--memory-limit", type=str, default="", help="Limit memory usage (e.g. '100G') in MMseqs2")

args = parser.parse_args()

fasta_path  = args.fasta_path
db_path = args.db_path
alignment_type = args.alignment_type
msa_out_dir = args.out_dir
n_cpu = args.n_cpu
n_iter = args.n_iter

if not args.max_hits:
    if alignment_type == "mgnify":
        max_hits = 501 
    elif alignment_type != "uniprot":
        max_hits = 10000
    else:
        max_hits = None
else:
    max_hits = args.max_hits

use_precomputed_msas = args.use_precomputed_msas

if alignment_type != "mmseqs2":
    # TODO: add mmseqs2 alignments
    if alignment_type in ["bfd_small", "uniref90", "mgnify", "uniprot"]:
        if not n_iter:
            n_iter = 1
        runner = jackhmmer.Jackhmmer(binary_path=which("jackhmmer"), database_path=db_path, n_cpu=n_cpu, n_iter=n_iter)
    elif alignment_type in ["uniref30", "bfd"]:
        if not n_iter:
            n_iter = 3
        runner = hhblits.HHBlits(binary_path=which("hhblits"), databases=[db_path], n_cpu=n_cpu, n_iter=n_iter)
        

    msa_format = "a3m" if alignment_type in ["uniref30", "bfd"] else "sto"
    msa_out_file = f"{alignment_type}_hits.{msa_format}"
    msa_out_path = os.path.join(msa_out_dir, msa_out_file)

    result = run_msa_tool(
            msa_runner=runner,
            input_fasta_path=fasta_path,
            msa_out_path=msa_out_path,
            msa_format=msa_format,
            use_precomputed_msas=use_precomputed_msas,
            max_sto_sequences=max_hits)

else:
    #sys.path.insert(0, args.colab_dir)
    #from colabfold.mmseqs import search

    max_accept = "--max-accept {max_hits}" if max_hits else ""
    mmseqs = Path(args.mmseqs) if args.mmseqs else which("mmseqs")

    Path(msa_out_dir).mkdir(parents=True, exist_ok=True)
    run_mmseqs(
        mmseqs,
        ["createdb", fasta_path,  os.path.join(msa_out_dir, "qdb"), "--shuffle", "0"],
    )

    mmseqs_search_monomer(
        mmseqs=mmseqs,
        dbbase=Path(db_path),
        base=Path(msa_out_dir),
        use_env=args.use_env,
        use_templates=0,
        max_accept=max_accept,
        threads=n_cpu,
        filter=1,
        db_load_mode=args.db_load_mode,
        uniref_db=args.db1,
        metagenomic_db=args.db2,
        split_memory_limit=args.memory_limit,
    )

    fasta_path.unlink()
    run_mmseqs(args.mmseqs, ["rmdb", args.base.joinpath("qdb")])
    run_mmseqs(args.mmseqs, ["rmdb", args.base.joinpath("qdb_h")])

