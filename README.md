# Efficiently run all-vs-all dimer predictions with AF_cache and mmseqs2

## Installing AF_cache

### Stand-alone installation

If you don't wish to use a Docker image, the installation and setup procedure is as follows:

1. [Install Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-macos-linux--wsl)

2. Set up conda environment, install dependencies:

```bash
# clone this repository
git clone https://github.com/clami66/AF_cache.git
cd AF_cache/

# install requirements with mamba/conda
mamba env create --file=environment.yaml
mamba activate AF_cache

# install python requirements
python -m pip install -r requirements.txt
```

## Nextflow PPI pipeline

1. Activate conda env:

```
conda activate AF_cache
export AF_CACHE=$(pwd)
```

2. Concatenate all fasta sequences into a single file. This is necessary only to run MMseqs2

```
$ cat fasta_seqs/*.fasta > all_seqs.fasta
```

3. Run the Nextflow workflow on a GPU node with e.g. 4 GPUs:

```
# -resume avoids re-running completed steps if the job crashed
$ nextflow AF_cache.nf -resume --fasta all.fasta --use_env --n_gpu 4 --proj_id berzelius-2025-6
```

To restrict the interactions to a list of `prot1_prot2` pairs, pass the absolute path to the list `multimers_list`:

```
$ head multimers_list 
YP00901869113_YP00901869113
YP00901869012_YP00901869012
...

$ nextflow AF_cache.nf -resume --fasta all.fasta --use_env --n_gpu 4 --proj_id berzelius-2025-6 --file_list $(realpath multimers_list)
```

## Running alignments and other steps separately

1. Run alignments with MMseqs2

GPU:
```
conda activate /proj/beyondfold/apps/.conda/envs/AF_cache
mmseqs_db=/proj/beyondfold/apps/colabfold_databases/gpu/
mmseqs_bin=/proj/beyondfold/apps/MMseqs2/build/bin/mmseqs

# 128 cpus for 4 GPUs on berzelius, 256 for 8 GPUs etc
$ python $AF_CACHE/run_msa_tool.py all.fasta mmseqs2 $mmseqs_db --out_dir align_outdir// --gpu --mmseqs $mmseqs_bin --n_cpu $n_cpu --use-env --n_cpu 128
```

CPU (not recommended):
```
conda activate /proj/beyondfold/apps/.conda/envs/AF_cache
mmseqs_db=/proj/beyondfold/apps/colabfold_databases/cpu/
mmseqs_bin=/proj/beyondfold/apps/MMseqs2/build/bin/mmseqs

$ python $AF_CACHE/run_msa_tool.py all.fasta mmseqs2 $mmseqs_db --out_dir align_outdir/ --mmseqs $mmseqs_bin --n_cpu $n_cpu --use-env --n_cpu 32
```


    Expected outputs:
    * align_outdir/alignments/ : folder including N .a3m files, one for each of the N fasta sequences in all_seqs.fasta


2. [optional] The alignments are written as a series of `.a3m` files in path/to/alignments/outdir/. These have to be converted to AlphaFold-like alignments, then we make AlphaFold parse the alignments and save them into pickle files. This is done in parallel for all a3m files at once:

```
$ python $AF_CACHE/prepare_alignments.py align_outdir/alignments/ AF_outdir/
$ ls fasta_seqs/*.fasta | parallel -j $n_cpu python $AF_CACHE/parse_features.py --flagfile $AF_CACHE/multimer_full_dbs_v3.flag --output_dir AF_outdir/ --fasta_paths={} --pickle_cache pickle_cache/ --alignments_only
```

    Expected outputs:

    * pickle_cache/ : one `.pkl` file per fasta sequence in the dataset

3. [optional] Generate all-vs-all fasta files, AF folder structures and package together jobs in multimer scripts:

```
$ python $AF_CACHE/format_alphafold_jobs.py fasta_seqs/ AF_outdir/ --pickle_dir pickle_cache --write_fastas --proj_id <berzelius-proj-id>
```

    Expected outputs:

    * AF_outdir/sbatch_scripts/ : all the packaged sbatch job scripts
        (e.g. /proj/beyondfold/users/x_clami/mmseqs_benchmark/egg_sperm_all/AF_models/sbatch_scripts/4500_0.sh)
    * AF_outdir/outdir/logs : all the slurm logs will go here
    * AF_outdir/<p1>_<p2> : all the dimers output folders. These will include files such as:
        * path/to/AF/outdir/<p1>_<p2>/<p1>_<p2>.fasta
        * path/to/AF/outdir/<p1>_<p2>/ranked_0.pdb
        * etc.

