#!/usr/bin/env nextflow

outputDir = '/proj/beyondfold/apps/alphafoldv2.3.1_pad/outputs'
params.conda_env = '/proj/beyondfold/apps/.conda/envs/AF_cache'
params.executor = 'slurm'
params.proj_id = 'berzelius-2024-439'
params.mmseqs_db = '/proj/common-datasets/LocalColabFold/'
params.mmseqs_bin = '/proj/beyondfold/apps/MMseqs2/build/bin/mmseqs'
params.af_dir = '/proj/beyondfold/apps/alphafoldv2.3.1_pad'
params.af_flagfile = '/proj/beyondfold/apps/alphafoldv2.3.1_pad/multimer_full_dbs_v3.flag'
params.use_env = false
params.n_gpu = 8
params.max_cpus = 64
params.file_list = ''
params.test = false


process split_fasta {
    executor = 'local'
    publishDir {outputDir}

    input:
    path fasta

    output:
    path "split_fasta/"

    script:
    """
    mkdir -p split_fasta/
    python ${params.af_dir}/split_fasta.py $fasta split_fasta/
    """
}

process ln_fasta {
    executor = 'local'
    publishDir {outputDir}

    input:
    path fasta

    output:
    path "fasta/*.fasta"

    script:
    """
    ln -s $fasta fasta
    """
}

process mmseqs_align {
    executor = "${params.executor}"
    clusterOptions "-A ${params.proj_id} -t 8:00:00 --gpus ${params.n_gpu}"
    publishDir {outputDir}
    
    input:
    path fasta
    
    output:
    path "alignments/"

    script:
    def use_env = params.use_env ? "--use-env" : ''
    def n_cpu = params.n_gpu > 0 ? 32*params.n_gpu : params.max_cpus
    """
    if ${params.test}; then
        mkdir -p alignments
        cp /proj/beyondfold/apps/alphafoldv2.3.1_pad/alignment_examples/*.a3m alignments/
    else
        python ${params.af_dir}/run_msa_tool.py $fasta mmseqs2 ${params.mmseqs_db} --out_dir ./ --gpu --mmseqs ${params.mmseqs_bin} --n_cpu $n_cpu $use_env
    fi
    """
}

process convert_alignments {
    executor = 'local'
    publishDir {outputDir}

    input:
    path alignments

    output:
    path "AF_data/"

    script:
    """
    python ${params.af_dir}/prepare_alignments.py $alignments AF_data/
    """
}


process parse_features {
    executor = 'local'
    publishDir {outputDir}, mode: 'copy'
    
    input:
    path fasta
    path af_data
    
    output:
    path "pickle_cache/*.gz"
    
    script:
    """
    mkdir -p pickle_cache
    python ${params.af_dir}/parse_features.py --flagfile ${params.af_flagfile} --output_dir $af_data --fasta_paths $fasta --pickle_cache pickle_cache/ --alignments_only
    """
}


process format_af_jobs {
    executor = 'local'
    publishDir {outputDir}

    input:
    path fasta
    path pickle_cache

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "AF_data_multimer/**.sh", emit: sh

    script:
    def file_list = params.file_list != '' ? "--file_list ${params.file_list}" : ''
    """
    python ${params.af_dir}/format_alphafold_jobs.py $fasta AF_data_multimer/ --conda_env ${params.conda_env} --pickle_dir ${outputDir}/pickle_cache --write_fastas --proj_id ${params.proj_id} $file_list
    """
}

process run_af_jobs {
    executor = 'slurm'
    clusterOptions "-A ${params.proj_id} -t 3-0 --gpus 1"
    publishDir {outputDir}
    
    input:
    path sbatch_script
    
    script:
    """
    sh $sbatch_script
    """
}

workflow {
    // align
    split_fasta_path = split_fasta(file(params.fasta))
    alignments_path = mmseqs_align(file(params.fasta))
    
    // convert
    af_data_path = convert_alignments(alignments_path)
    pickle_cache = parse_features(ln_fasta(split_fasta_path).flatten(), af_data_path).collect().flatten().take(1)
    
    // AF
    sbatch_scripts = format_af_jobs(split_fasta_path, pickle_cache).sh.collect().flatten()
    run_af_jobs(sbatch_scripts)
}
