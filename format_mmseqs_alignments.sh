#!/bin/bash
set -euxo pipefail

AFPATH=/proj/beyondfold/apps/alphafoldv2.3.1_pad/
DBDIR=/proj/beyondfold/apps/colabfold_databases/
FLAGFILE=/proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/multimer_all_vs_all.flag

a3m_path=$(realpath $1)
fastapath=$2
outdir=$3
picklecache=$4

mkdir -p $picklecache

#python $AFPATH/mmseqs_2_uniprot_parallel.py $a3m_path $outdir $DBDIR --n_cpu 16 --shortname
ls $fastapath/*.fasta | parallel -j 16 python $AFPATH/parse_features.py --flagfile $FLAGFILE --output_dir $outdir --fasta_paths={} --pickle_cache $picklecache --alignments_only
