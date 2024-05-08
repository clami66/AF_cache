#!/bin/bash
set -euxo pipefail

AFPATH=/proj/beyondfold/apps/alphafoldv2.3.1_pad/
DBDIR=/proj/beyondfold/apps/colabfold_databases/
FLAGFILE=/proj/beyondfold/apps/alphafoldv2.3.1_pad/multimer_full_dbs_v3.flag

a3m_path=$(realpath $1)
fastapath=$2
outdir=$3
picklecache=$4

mkdir -p $picklecache

python $AFPATH/mmseqs_2_uniprot_parallel.py $a3m_path $outdir $DBDIR --n_cpu 16
ls $fastapath/*.fasta | parallel -j 16 python $AFPATH/run_alphafold.py --flagfile $FLAGFILE --output_dir $outdir --fasta_paths={} --pickle_cache $picklecache --alignments_only
