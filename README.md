# AF_cache: fast inference of AlphaFold2 / AlphaFold3 predictions in large-scale studies

## Setup

1. [Install nextflow](https://www.nextflow.io/docs/latest/install.html)

    ```bash
    curl -s https://get.nextflow.io | bash
    mv nextflow ~/bin/
    ```

2. Clone this repository and export the install path as follows:
   ```bash
   git clone git@github.com:clami66/AF_cache.git
   cd AF_cache/

   # export installation path to bash profile, change accordingly if using e.g. zsh
   echo "export AF_CACHE=$(pwd)" >> ~/.bashrc
   source ~/.bashrc
   ```

3. Run the pipeline for the first time. This will automatically download and setup all the necessary DBs/tools. This could take a few hours the first time you run the pipeline.
    ```
    # add --test to skip the alignment step
    nextflow AF_cache.nf --fasta test_data/fasta/all.fasta -resume
    ```

* **The pipeline uses Docker containers to automatically get all the requirements. Alternatively, apptainer or conda can also be used. See below to configure this behavior.**

* **The pipeline automatically installs ColabFold MSA DBs and AlphaFold2 PDB template DBs. If these are already present on the system, this step can be skipped. See below to configure this behavior.**

### Pipeline inputs

The input to a pipeline is a single `.fasta` file containing all the sequences for a large-scale experiment:

```
cat fasta_seqs/*.fasta > all_seqs.fasta
```

The workflow will take this single fasta file as input:

```
# -resume avoids re-running completed steps if the job crashed
nextflow AF_cache.nf --fasta all_seqs.fasta -resume
```

### Subsetting pairs, running multimers (trimers, etc.)

By default, the pipeline produces predictions for all-vs-all pairs of sequences in the input fasta.

To restrict the interactions to a list of `prot1 prot2` pairs, use the parameter `--pair_list`. Even though the pipeline is thought for dimer interactions, more than two partners can be specified:
    
```
$ head multimers_list 
YP00901869113 YP00901869113
YP00901869012 YP00901869012 YP00901869113
...

# will generate predictions for one homodimer and one heterotrimer:
nextflow AF_cache.nf --fasta all.fasta --pair_list multimers_list
```

## Other configuration options

The pipeline's behavior can be customized depending on what is available on the host system. Most configuration is done within `nextflow.config`:

<details>
<summary>Choosing between Docker, apptainer/singularity and conda</summary>

The pipeline uses Docker containers to automatically get all the requirements. Alternatively, apptainer or conda can also be used by running different profiles:

Apptainer/singularity:
```
nextflow AF_cache.nf --fasta test_data/fasta/all.fasta -profile apptainer
```
Conda/Mamba:
```
nextflow AF_cache.nf --fasta test_data/fasta/all.fasta -profile conda
```

</details>   
<details>
<summary>Running AlphaFold3</summary>
AlphaFold3 can be run by simply adding the `--af3` flag:

```
nextflow AF_cache.nf --fasta test_data/fasta/all.fasta --af3
```
This will automatically install the necessary environment, according to the docker/apptainer/conda preferences described above.

**Notice: the AF3 docker container is not maintained by us.**
</details>

<details>
<summary>Skipping templates</summary>
The template step in the pipeline can be enabled/skipped, either by permanently setting `skip_template` inside `nextflow.config`:

```
skip_templates = true
```

By default, templates are always skipped. This will also avoid downloading the template DBs the first time the pipeline is run.

Like all `params` options, this behavior can be changed at runtime:
```
nextflow AF_cache.nf --fasta test_data/fasta/all.fasta --skip_templates=false
```
</details>

<details>
<summary>ColabFold DBs are already on the system</summary>
If the ColabFold DBs have been downloaded through the original ColabFold setup script, these can be reused so that the pipeline doesn't download an extra copy. This can be done by pointing the `database` directory inside `nextflow.config` to the right location:

```
mmseqs_db = "/path/to/ColabFold/DB"
```
The directory should contain the files `DOWNLOADS_READY`, `UNIREF30_READY`, `COLABDB_READY`. The pipeline will look for these files to skip the download step.
</details>

<details>
<summary>AlphaFold2 template DBs are already on the system</summary>
If template DBs (`pdb_mmcif`, `pdb_seqres`) are already on the system, these can be used to avoid downloading an extra copy. 
    
1. Set the correct paths inside `nextflow.config`
```
template_mmcif_dir = "/path/to/pdb_mmcif/mmcif_files"
obsolete_pdbs_path = "/path/to/pdb_mmcif/obsolete.dat"
pdb_seqres_database_path = "/path/to/pdb_seqres/pdb_seqres.txt"
```
2. Make sure you add an empty file called `PDB_MMCIF_READY` inside the directory of the ColabFold DBs (`mmseqs_db = "/path/to/ColabFold/DB"`). This will avoid triggering a new download of the DBs.
   
</details>

<details>
<summary>AlphaFold2/AlphaFold3 configuration</summary>

AlphaFold2 and/or AlphaFold3 parameters should be on the system, according to the respective installation instructions. Simply point the pipeline to the parameter file locations inside `nextflow.config`:

```
af2_data_dir = '/path/to/alphafold2_data/'
af3_model_dir = '/path/to/af3_model_parameters/'
```

Other behaviors for AF2 and AF3 (number of recycles, number of seeds etc.) should be set inside the flagfiles provided inside `flags/af2.flag` and `flags/af3.flag`. For example, to change the number of recycles inside AF2 and use two NN models, `flags/af2.flag` might look as follows:

```
--max_recycles=3
--models_to_use=model_1_multimer_v3,models_2_multimer_v3
```

</details>
<details>

<summary>Job scheduling and resource management</summary>

Depending whether the pipeline runs on an HPC sytem or locally, some parameters can be varied to send jobs to different schedulers or to run them on local CPU/GPU resourcers.

For example: on a SLURM-based system, one could send the alignment job to a node with 8 GPUs and all AF inference jobs to single-GPU nodes. Other lighter tasks (e.g. parsing features, copying files) can be sent to CPU-only nodes, or run locally (on the front node). That would be accomplished with the following settings in `nextflow.config`:

```
withName:'mmseqs_align' {
    executor = 'slurm'
    clusterOptions = '--account xxx-yyy-zzz --gpus 8 --time 12:00:00'
}

withName:'run_af2_jobs|run_af3_jobs' {
    executor = 'slurm'
    clusterOptions = '--account xxx-yyy-zzz --gpus 1 --time 12:00:00'
}

withName:'ln_fasta|split_fasta|collect_pickles|collect_jsons' {
    executor = 'local'
}

withName:'convert_alignments_af2|convert_alignments_af2|parse_features_af2|parse_features_af3|format_jobs_af2|format_jobs_af3' {
    executor = 'local'
}
```

If all tasks are running on a local machine, all executors may be set to `local`, then edit the executor to make sure that the job queue size is not larger than the number of GPUs available on said machine. For example, if four GPUs are on a local machine:

```
executor{
    name = "local"
    queueSize = 4
    cpus = 32
}
```

Heavier CPU tasks to process alignments may be sent to a CPU node, e.g. through SLURM:

```
...

withName:'convert_alignments_af2|convert_alignments_af2|parse_features_af2|parse_features_af3|format_jobs_af2|format_jobs_af3' {
    executor = 'slurm'
    clusterOptions = '--account xxx-yyy-zzz -N 1 -n 32 --time 2:00:00'
}
```

Consult the Nextlow docs for more information about setting up different executors/schedulers [here](https://www.nextflow.io/docs/latest/reference/config.html#executor).

</details>

## Running pipeline steps separately (deprecated)

1. Run alignments with MMseqs2

GPU:
```
conda activate AF_cache
mmseqs_db=/path/to/mmseqs_db
mmseqs_bin=/path/to/mmseqs

python $AF_CACHE/run_msa_tool.py all.fasta mmseqs2 $mmseqs_db --out_dir align_outdir/ --gpu --mmseqs $mmseqs_bin --n_cpu $n_cpu --use-env --n_cpu 128
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
$ python $AF_CACHE/format_alphafold_jobs.py fasta_seqs/ AF_outdir/ --pickle_dir pickle_cache --write_fastas --proj_id <proj-id>
```

    Expected outputs:

    * AF_outdir/sbatch_scripts/ : all the packaged sbatch job scripts
        (e.g. /proj/beyondfold/users/x_clami/mmseqs_benchmark/egg_sperm_all/AF_models/sbatch_scripts/4500_0.sh)
    * AF_outdir/outdir/logs : all the slurm logs will go here
    * AF_outdir/<p1>_<p2> : all the dimers output folders. These will include files such as:
        * path/to/AF/outdir/<p1>_<p2>/<p1>_<p2>.fasta
        * path/to/AF/outdir/<p1>_<p2>/ranked_0.pdb
        * etc.

