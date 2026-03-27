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
# export the following to make sure that the default pipeline config settings are correct
export AF_CACHE=$(pwd)

# install requirements with mamba/conda
mamba env create --file=environment.yaml
mamba activate AF_cache

# install python requirements
python -m pip install -r requirements.txt
```
3. Download and setup MMseqs2 (more detailed instructions [here](https://github.com/soedinglab/mmseqs2))
```bash
# GPU version:
wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz; tar xvfz mmseqs-linux-gpu.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

# Non-GPU version:
# wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
```

5. Download and setup ColabFold DBs (instructions taken [here](https://colabfold.mmseqs.com/))

```bash
wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
chmod +x setup_databases.sh
# If you want to use GPU acceleration during the searches:
GPU=1 ./setup_databases.sh database/
# If you DO NOT want to use GPU acceleration during the searches:
# ./setup_databases.sh database/
```

## Nextflow PPI pipeline (recommended)

### Configuring the pipeline

Most configuration is done in `nextflow.config`. Here, one can set up the installation paths of AF_cache, the ColabFold DBs, MMseqs2.

#### Installation paths

If the installation instructions were followed exactly, there is no need to change these:

```
conda_env = 'AF_cache'
# necessary to export AF_CACHE install path
af_cache_dir = '$AF_CACHE'
mmseqs_db = '$AF_CACHE/database/'
mmseqs_bin = '$AF_CACHE/mmseqs/bin/mmseqs'
```

If the DBs and MMseqs2 were installed in some other location, the parameters need to be adjusted accordingly.

#### AlphaFold parameters

Other parameters can be adjusted to change the behavior of AlphaFold2, or to point to an existing installation of AlphaFold3 if the user would like to run `AF3_cache.nf` instead.

```
// af2 parameters
af_flagfile = '$AF_CACHE/flags/multimer.flag'    // flagfile with default AF2.3 inference parameters
db_flagfile = '$AF_CACHE/flags/databases.flag'   // flagfile with DB paths for AF2.3. Might have to be adjusted to point to DB location

// af3 parameters
af3_conda_env = 'your_AF3_conda_environment'
af3_dir = '/path/to/your/af3/installation/'
```

**NB:** the flagfiles (`af_flagfile`, `db_flagfile`, etc.) are a convenient place to set all the necessary flags to run AlphaFold. In this repo, we have a set of predefined flagfiles (inside `flags/`). These need to be adjusted, for example, so that AF can find the model parameters and the ColabFold databases. It is important to make sure that the information in the flagfiles are correct.

#### Scheduling and resource management

Depending whether the pipeline runs on an HPC sytem or locally, some parameters can be varied to send jobs to different schedulers or to run them on local CPU/GPU resourcers.

For example: on a SLURM-based system, one could send the alignment job to a node with 8 GPUs and all AF inference jobs to single-GPU nodes. Other lighter tasks (e.g. parsing features, copying files) can be sent to CPU-only nodes, or run locally (on the front node). That would be accomplished with the following settings in `nextflow.config`:

```
mmseqs_executor = 'slurm'
mmseqs_executor_flags = '--account your-account-ID --gpus 8 --time 12:00:00'

af_executor = 'slurm'
af_executor_flags = '--account your-account-ID --gpus 1 --time 12:00:00'

other_executor = 'local'
other_executor_flags = ''
```

If all tasks are running on a local machine, one can set all executors to `local`, then edit the executor to make sure that the job queue size is not larger than the number of GPUs available on said machine. For example, if four GPUs are on a local machine:

```
executor{
    name = "local"
    queueSize = 4
    cpus = 32
}
```

Consult the Nextlow docs for more information about setting up different executors/schedulers [here](https://www.nextflow.io/docs/latest/reference/config.html#executor).

### Running the pipeline

1. Activate conda env:

```
conda activate AF_cache
```

2. Concatenate all fasta sequences into a single file. This is necessary only to run MMseqs2

```
$ cat fasta_seqs/*.fasta > all_seqs.fasta
```

3. Run the Nextflow workflow on a GPU node with e.g. 4 GPUs:

```
# -resume avoids re-running completed steps if the job crashed
$ nextflow AF_cache.nf -resume --fasta all.fasta --use_env --n_gpu 4 --proj_id slurm-proj-id
```

To restrict the interactions to a list of `prot1 prot2` pairs, pass the absolute path to the list `multimers_list`:

```
$ head multimers_list 
YP00901869113 YP00901869113
YP00901869012 YP00901869012
...

$ nextflow AF_cache.nf -resume --fasta all.fasta --use_env --n_gpu 4 --proj_id slurm-proj-id --file_list $(realpath multimers_list)
```

## Running alignments and other steps separately (not recommended)

1. Run alignments with MMseqs2

GPU:
```
conda activate AF_cache
mmseqs_db=/path/to/mmseqs_db
mmseqs_bin=/path/to/mmseqs

# 128 cpus for 4 GPUs on berzelius, 256 for 8 GPUs etc
$ python $AF_CACHE/run_msa_tool.py all.fasta mmseqs2 $mmseqs_db --out_dir align_outdir/ --gpu --mmseqs $mmseqs_bin --n_cpu $n_cpu --use-env --n_cpu 128
```

CPU (not recommended):
```
python $AF_CACHE/run_msa_tool.py all.fasta mmseqs2 $mmseqs_db --out_dir align_outdir/ --mmseqs $mmseqs_bin --n_cpu $n_cpu --use-env --n_cpu 32
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

