# Efficiently run all-vs-all dimer predictions with AF_cache and mmseqs2

## Installing AF_cache

1. [Install nextflow](https://www.nextflow.io/docs/latest/install.html)

    ```bash
    curl -s https://get.nextflow.io | bash
    mv nextflow ~/bin/
    ```

3. Clone this repository and export the install path as follows:
   ```bash
   git clone https://github.com/clami66/AF_unmasked.git
   cd AF_cache/
   export AF_CACHE=$(pwd)
   ```

4. Download and setup ColabFold DBs (instructions taken [here](https://colabfold.mmseqs.com/))

    ```bash
    wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh
    chmod +x setup_databases.sh
    # If you want to use GPU acceleration during the searches:
    GPU=1 ./setup_databases.sh database/
    # If you DO NOT want to use GPU acceleration during the searches:
    # ./setup_databases.sh database/
    ```

5. Install the environments necessary to run AlphaFold and MMseqs2. This can be done in multiple ways:
   <details>
   <summary>Using Docker</summary>

   If you wish to use our Docker container, just open `nextflow.config` in the main repository folder and enable `docker` while disabling `apptainer` and `conda`:
   ```
   docker {
        enabled = true
        runOptions = "--gpus all"
    }
   
   apptainer {
        enabled = false
        runOptions = "--nv"
    }

    conda {
        enabled = false
        useMamba = true
        cacheDir = "$AF_CACHE/conda"
    }
   ```

   Then, run the pipeline for the first time and the correct image will be pulled by nextflow.
   </details>
   <details>
   <summary>Using Apptainer</summary>

   If you wish to use Apptainer, just open `nextflow.config` in the main repository folder and enable `apptainer` while disabling `docker` and `conda`:
   ```
   docker {
        enabled = false
        runOptions = "--gpus all"
    }
   
   apptainer {
        enabled = true
        runOptions = "--nv"
    }

    conda {
        enabled = false
        useMamba = true
        cacheDir = "$AF_CACHE/conda"
    }
   ```

   Then, run the pipeline for the first time and the correct image will be pulled by nextflow. You can specify where the `.sif` container image will be downloaded by setting the `NXF_APPTAINER_CACHEDIR` ennvar as desired.
   </details>
   <details>
   <summary>Using conda/mamba</summary>

   If you wish to use conda or mamba, a few extra steps are needed. First, you need to install MMseqs2 (for GPU or CPU) as follows:
   ```bash
   # GPU version:
   wget https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz; tar xvfz mmseqs-linux-gpu.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

   # Non-GPU version:
   # wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
   ```

   Then, add the path to the `mmseqs` binary inside `nextflow.config`. If you install MMseqs2 inside the main repo directory, it should look as follows:

   ```
   mmseqs_bin = '$AF_CACHE/mmseqs/bin/mmseqs'
   ```

   Lastly, enable conda inside `nextflow.config`:
   ```
   docker {
        enabled = false
        runOptions = "--gpus all"
    }
   
   apptainer {
        enabled = false
        runOptions = "--nv"
    }

    conda {
        enabled = true
        useMamba = true
        cacheDir = "$AF_CACHE/conda"
    }
   ```
   Then, run the pipeline for the first time and the conda environment will be automatically set up.
   You can also specify where the environment should be installed by setting `cacheDir`.
   </details>   



6. Test the pipeline

    ```
    # --test skips the alignment step
    nextflow AF_cache.nf --fasta test_data/fasta/all.fasta --test -resume
    ```

### Configuring the pipeline

Most configuration is done within `nextflow.config`. Here, one can set up the installation paths of AF_cache, the ColabFold DBs, MMseqs2.

#### Installation paths

If the installation instructions were followed exactly, there is no need to change these:

```
# necessary to export AF_CACHE install path
af_cache_dir = '$AF_CACHE'
mmseqs_db = '$AF_CACHE/database/'
mmseqs_bin = '$AF_CACHE/mmseqs/bin/mmseqs'
```

If the DBs and MMseqs2 were installed in some other location, the parameters need to be adjusted accordingly.

#### AlphaFold parameters

Other parameters can be adjusted to change the behavior of AlphaFold2, or to point to an existing installation of AlphaFold2/AlphaFold3. Notice that the original AF2 template databases are still needed (`template_mmcif_dir`, `obsolete_pdbs_path`, `pdb_seqres_database_path`).

```
// af2 parameters
af2_flagfile = '$AF_CACHE/flags/multimer.flag'
af2_data_dir = '/path/to/af2/parameters'
template_mmcif_dir = '.../pdb_mmcif/mmcif_files'
obsolete_pdbs_path = '.../pdb_mmcif/obsolete.dat'
pdb_seqres_database_path = '.../pdb_seqres.txt'

// af3 parameters
af3_dir = '/path/to/your/af3/installation/'
...
```

**NB:** the flagfiles (`af_flagfile`, `db_flagfile`, etc.) are a convenient place to set all the necessary flags to run AlphaFold. In this repo, we have a set of predefined flagfiles (inside `flags/`).

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

1. Concatenate all fasta sequences into a single file:

    ```
    $ cat fasta_seqs/*.fasta > all_seqs.fasta
    ```

3. Launch the Nextflow workflow:

    ```
    # -resume avoids re-running completed steps if the job crashed
    $ nextflow AF_cache.nf -resume --fasta all.fasta
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

