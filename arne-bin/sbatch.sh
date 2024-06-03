#!/bin/bash
#SBATCH -A berzelius-2023-295
#SBATCH --gpus 1
#SBATCH -t 3-00:00:00
#SBATCH -C thin
#SBATCH -o logs/afrun.log

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='6.0'

module load Anaconda/2021.05-nsc1

conda activate /proj/beyondfold/apps/.conda/envs/af_server 



for id in P57053_P04487 P62937_P04487 Q15056_P10225
do 
    cache=cache_test/
    outdir=multimer_test/
    python /proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all.flag --output_dir $outdir --fasta_paths $outdir/$id/$id.fasta --pickle_cache $cache  2> logs/$id.default.log


    cache=cache_test/identical/
    outdir2=multimer_test_identical/
    python /proj/berzelius-2021-29/users/x_arnel/herpes/herpes-ppi/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all.flag --output_dir $outdir2 --fasta_paths $outdir/$id/$id.fasta --pickle_cache $cache 2> logs/$id.includeidentical.log
done

# python /proj/beyondfold/apps/alphafoldv2.3.1_cache/run_alphafold.py --flagfile /proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/multimer_all_vs_all.flag --output_dir /proj/berzelius-2021-29/users/x_arnel/histone/multimer --fasta_paths /proj/berzelius-2021-29/users/x_arnel/histone/multimer/Hs-Samp1a_Histone2b/Hs-Samp1a_Histone2b.fasta --pickle_cache cache/ --pad_to_size 800,20000 
