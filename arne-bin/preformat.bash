#!/bin/bash -x

module load Mambaforge/23.3.1-1-hpc1-bdist
conda activate /proj/beyondfold/apps/.conda/envs/af_server

#bin/runmmseq.bash all.fasta
#python /proj/beyondfold/apps/alphafoldv2.3.1_pad/run_msa_tool.py seq/all.fasta mmseqs2 /proj/beyondfold/apps/colabfold_databases/ --out_dir msa/ --mmseqs /proj/beyondfold/apps/.conda/envs/colabfold/bin/mmseqs

rm -rf cache*/
rm -rf multimer*/

bin/format_mmseqs_alignments.sh msa/ seq/ multimer/ cache/ /proj/beyondfold/apps/colabfold_databases/
python bin/format_alphafold_jobs_AE.py seq/ multimer/  --pickle_dir  cache/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --max_job_size 250 100 25 20 10 5 1
python bin/format_alphafold_jobs_AE.py seq/ multimer25/  --pickle_dir  cache/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25.flag --max_job_size 100 50 10 10 10 5 1

bin/format_mmseqs_alignments_duplicate.sh msa/ seq/ multimer_family_dupl/ cache_dupl/family/ family_duplicate/
python bin/format_alphafold_jobs_AE.py seq/ multimer_family_dupl/  --pickle_dir  cache_dupl/family/  --write_fastas  --proj_id berzelius-2023-295  --file_list pairs.txt  --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all.flag --max_job_size 250 100 25 20 10 5 1
python bin/format_alphafold_jobs_AE.py seq/ multimer_family_dupl25/  --pickle_dir  cache_dupl/family/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25.flag --max_job_size 100 50 10 10 10 5 1


rsync -arv multimer_family_dupl/ multimer_family_id/
#bin/format_mmseqs_alignments.sh msa/ seq/ multimer_family_id/ cache_id/family/ family_duplicate/
python bin/format_alphafold_jobs_AE.py seq/ multimer_family_id/  --pickle_dir  cache_dupl/family/identical/  --write_fastas  --proj_id berzelius-2023-295  --file_list pairs.txt  --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_includeidentical.flag --max_job_size 250 100 25 20 10 5 1
python bin/format_alphafold_jobs_AE.py seq/ multimer_family_id25/  --pickle_dir  cache_dupl/family/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25_identical.flag --max_job_size 100 50 10 10 10 5 1




bin/format_mmseqs_alignments_duplicate.sh msa/ seq/ multimer_subfamily_dupl/ cache_dupl/subfamily/ subfamily_duplicate/
python bin/format_alphafold_jobs_AE.py seq/ multimer_subfamily_dupl/  --pickle_dir  cache_dupl/subfamily/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --max_job_size 250 100 25 20 10 5 1
#rsync -arv multimer_subfamily_dupl/ multimer_subfamily_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_subfamily_id/  --pickle_dir  cache_dupl/subfamily//identical  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_includeidentical.flag --max_job_size 250 100 25 20 10 5 1

bin/format_mmseqs_alignments_duplicate.sh msa/ seq/ multimer_genus_dupl/ cache_dupl/genus/ genus_duplicate/
python bin/format_alphafold_jobs_AE.py seq/ multimer_genus_dupl/  --pickle_dir  cache_dupl/genus/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --max_job_size 250 100 25 20 10 5 1
#rsync -arv multimer_genus_dupl/ multimer_genus_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_genus_id/  --pickle_dir  cache_dupl/genus/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_includeidentical.flag --max_job_size 250 100 25 20 10 5 1

bin/format_mmseqs_alignments_duplicate.sh msa/ seq/ multimer_order_dupl/ cache_dupl/order/ order_duplicate/
python bin/format_alphafold_jobs_AE.py seq/ multimer_order_dupl/  --pickle_dir  cache_dupl/order/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt  --max_job_size 250 100 25 20 10 5 1
#rsync -arv multimer_order_dupl/ multimer_order_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_order_id/  --pickle_dir  cache_dupl/order/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_includeidentical.flag --max_job_size 250 100 25 20 10 5 1


bin/format_mmseqs_alignments_duplicate.sh msa/ seq/ multimer_sprot_dupl/ cache_sprot/ SwissProt_dupl/ 
python bin/format_alphafold_jobs_AE.py seq/ multimer_sprot_dupl/  --pickle_dir  cache_sprot/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt  --max_job_size 250 100 25 20 10 5 1
#rsync -arv multimer_sprot_dupl/ multimer_sprot_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_sprot_id/  --pickle_dir  cache_sprot/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_includeidentical.flag  --max_job_size 250 100 25 20 10 5 1

python bin/format_alphafold_jobs_AE.py seq/ multimer_sprot_dupl25/  --pickle_dir  cache_sprot/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt  --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25.flag  --max_job_size 40 20 10 10 10 5 1
#rsync -arv multimer_sprot_dupl/ multimer_sprot_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_sprot_id25/  --pickle_dir  cache_sprot/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25_identical.flag  --max_job_size 40 20 10 10 10 5 1

bin/format_mmseqs_alignments_duplicate.sh msa/ seq/ multimer_trembl_dupl/ cache_trembl/ TremblSprot/ 
python bin/format_alphafold_jobs_AE.py seq/ multimer_trembl_dupl/  --pickle_dir  cache_trembl/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt  --max_job_size 250 100 25 20 10 5 1
#rsync -arv multimer_sprot_dupl/ multimer_sprot_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_trembl_id/  --pickle_dir  cache_trembl/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_includeidentical.flag  --max_job_size 250 50 20 15 10 5 1

python bin/format_alphafold_jobs_AE.py seq/ multimer_trembl_dupl25/  --pickle_dir  cache_trembl/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt  --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25.flag  --max_job_size 40 10 8 5 3 2 1
#rsync -arv multimer_sprot_dupl/ multimer_sprot_id/
python bin/format_alphafold_jobs_AE.py seq/ multimer_trembl_id25/  --pickle_dir  cache_trembl/identical/  --write_fastas  --proj_id berzelius-2023-295   --file_list pairs.txt --flagfile /proj/berzelius-2021-29/users/x_arnel/herpes/bin/multimer_all_vs_all_25_identical.flag  --max_job_size 40 20 10 10 10 5 1



# Generating analysis data
find ./ -size 0 -name "*.csv" -exec rm {} \; # remove empty files
nohup ls -d  multimer*/*_*/  | parallel -j 16 ../bin/extractall.bash {} &
find ./ -size 0 -name "*.csv" -exec rm {} \; # remove empty files

#for i in multimer* ; do for j in pdockq pdockq_v21 ptm  ; do cat $i/*_*/$j.csv | grep -v sbatch > herpes-ppi/data/$i-$j.csv ; done ; done