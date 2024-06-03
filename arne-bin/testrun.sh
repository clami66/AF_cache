#!/bin/bash
#SBATCH -A berzelius-2023-295
#SBATCH --gpus 1
#SBATCH -t 0-01:00:00
#SBATCH -C thin
#SBATCH -o /proj/berzelius-2021-29/users/x_arnel/herpes//logs/tesrun.log

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='6.0'

module load Anaconda/2021.05-nsc1

conda activate /proj/beyondfold/apps/.conda/envs/af_server 
rm -f /proj/berzelius-2021-29/users/x_arnel/herpes/multimer_test/P57053_P04487/*.pkl 
rm -f /proj/berzelius-2021-29/users/x_arnel/herpes/multimer_test/P57053_P04487/*.pdb
rm -f /proj/berzelius-2021-29/users/x_arnel/herpes/multimer_test/P57053_P04487/*.csv 
rm -f /proj/berzelius-2021-29/users/x_arnel/herpes/multimer_test/P57053_P04487/*.json
python /proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_featpickle.flag --output_dir /proj/berzelius-2021-29/users/x_arnel/herpes/multimer_test --fasta_paths /proj/berzelius-2021-29/users/x_arnel/herpes/multimer_test/P57053_P04487/P57053_P04487.fasta  --pickle_cache cache/ --pad_to_size 400,20000 
 
 python /proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_featpickle.flag --output_dir /proj/berzelius-2021-29/users/x_arnel/herpes/multimer --fasta_paths seq/P04487.fasta --pickle_cache cache/ --alignments_only 