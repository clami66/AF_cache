# Efficiently run all-vs-all dimer predictions with AF_pad and mmseqs2

1. Activate conda env:

```
module load Mambaforge/23.3.1-1-hpc1-bdist
conda activate /proj/beyondfold/apps/.conda/envs/af_server

# alternately you can just launch python scripts from outside the environment with:
# conda run -p /proj/beyondfold/apps/.conda/envs/af_server python ...
```

2. Concatenate all fasta sequences into a single file. This is necessary only to run MMseqs2

```
cat fasta_seqs/*.fasta > all_seqs.fasta
```

2. Run alignments with MMseqs2

```
mmseqs_db=/proj/beyondfold/apps/colabfold_databases/
mmseqs_bin=/proj/beyondfold/apps/.conda/envs/colabfold/bin/mmseqs
python /proj/beyondfold/apps/alphafoldv2.3.1_pad/run_msa_tool.py all_seqs.fasta mmseqs2 $mmseqs_db --out_dir path/to/alignments/outdir --mmseqs $mmseqs_bin
```

    Expected outputs:
    * path/to/alignments/outdir : folder including N .a3m files, one for each of the N fasta sequences in all_seqs.fasta

3. The alignments are written as a series of `.a3m` files in path/to/alignments/outdir/. These have to be converted to AlphaFold-like alignments, then we make AlphaFold parse the alignments and save them into pickle files. This is done in parallel for all a3m files at once:

```
# first parameter: the alignment output dir from the previous step
# second parameter: the original folder containing all the fasta files
# third parameter: the main path where all AlphaFold predictions will be generated
# fourth parameter: a directory that will contain all of the feature `.pkl` files (one for each sequence in fasta dir)
/proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/format_mmseqs_alignments.sh path/to/alignments/outdir/ fasta_seqs/ path/to/AF/outdir path/to/pickle/cache
```

    Expected outputs:

    * path/to/pickle_cache/ : one `.pkl` file per fasta sequence in the dataset

3. Generate all-vs-all fasta files, AF folder structures and package together jobs in multimer scripts:

```
python /proj/beyondfold/users/x_clami/mmseqs_benchmark/scripts/format_alphafold_jobs.py fasta_seqs/ path/to/AF/outdir --pickle_dir path/to/pickle/cache --write_fastas --proj_id <berzelius-proj-id>
```

    Expected outputs:

    * path/to/AF/outdir/sbatch_scripts/ : all the packaged sbatch job scripts
        (e.g. /proj/beyondfold/users/x_clami/mmseqs_benchmark/egg_sperm_all/AF_models/sbatch_scripts/4500_0.sh)
    * path/to/AF/outdir/logs : all the slurm logs will go here
    * path/to/AF/outdir/<p1>_<p2> : all the dimers output folders. These will include files such as:
        * path/to/AF/outdir/<p1>_<p2>/<p1>_<p2>.fasta
        * path/to/AF/outdir/<p1>_<p2>/ranked_0.pdb
        * etc.

