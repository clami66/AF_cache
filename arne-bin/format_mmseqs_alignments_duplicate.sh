#!/bin/bash -x
set -euxo pipefail

#SCRIPTS=/proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/
SCRIPTS=/proj/berzelius-2021-29/users/x_arnel/herpes/bin/
#AFPATH=/proj/beyondfold/apps/alphafoldv2.3.1_pad/
AFPATH=/proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/
#DBDIR=/proj/beyondfold/apps/colabfold_databases/
DBDIR=$5
FLAGFILE=$SCRIPTS/multimer_all_vs_all.flag

a3m_path=$(realpath $1)
fastapath=$2
outdir=$3
picklecache=$4

mkdir -p $picklecache

python $SCRIPTS/mmseqs_2_uniprot_parallel_test.py $a3m_path $outdir $DBDIR --n_cpu 16 --duplicate


ls $fastapath/*.fasta | parallel -j 16 python $AFPATH/run_alphafold.py --flagfile $FLAGFILE --output_dir $outdir --fasta_paths={} --pickle_cache $picklecache --alignments_only
picklecache=$4/identical
mkdir -p $picklecache
FLAGFILE=$SCRIPTS/multimer_all_vs_all_includeidentical.flag
ls $fastapath/*.fasta | parallel -j 16 python $AFPATH/run_alphafold.py --flagfile $FLAGFILE --output_dir $outdir --fasta_paths={} --pickle_cache $picklecache --alignments_only
