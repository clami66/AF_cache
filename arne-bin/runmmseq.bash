#!/bin/bash -x
conda activate /proj/beyondfold/apps/.conda/envs/af_server

mmseqs_db=/proj/beyondfold/apps/colabfold_databases/
mmseqs_bin=/proj/beyondfold/apps/.conda/envs/colabfold/bin/mmseqs
#msadir=msa/
msadir=$2          
python /proj/beyondfold/apps/alphafoldv2.3.1_pad/run_msa_tool.py $1 mmseqs2 $mmseqs_db --out_dir $msadir --mmseqs $mmseqs_bin
