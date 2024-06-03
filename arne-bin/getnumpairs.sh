#!/bin/bash -x
#SBATCH -A berzelius-2023-295
#SBATCH --gpus 1
#SBATCH -t 3-00:00:00
#SBATCH -C thin
#SBATCH -o logs/afrun.log

#module load Anaconda/2021.05-nsc1

#conda activate /proj/beyondfold/apps/.conda/envs/af_server % This is the conda environment that has alphafoldv2.3.1_cache installed

AFPATH=/proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/run_alphafold.py

fasta=$1
cache=$2
outdir=`mktemp -d`
#id=$outdir/$id/$id.fasta
#flag=$2
flag=/proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all.flag
#flag=$3

python $AFPATH --flagfile $flag --output_dir $outdir --fasta_paths $fasta --pickle_cache $cache --alignments_only 2> $outdir/err.log
echo -n $fasta"," | sed "s/\//,/g" > $fasta.numpairs.csv
grep PAIRING $outdir/err.log | grep "Chain 0" | gawk '{print $9}' >> $fasta.numpairs.csv

rm -rf $outdir

