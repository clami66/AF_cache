#!/usr/bin/env nextflow

outputDir = '/proj/beyondfold/apps/alphafoldv2.3.1_pad/outputs'
params.conda_env = '/proj/beyondfold/apps/.conda/envs/af3'
params.executor = 'slurm'
params.proj_id = 'berzelius-2026-12'
params.mmseqs_db = '/proj/common-datasets/LocalColabFold/'
params.mmseqs_bin = '/proj/beyondfold/apps/MMseqs2/build/bin/mmseqs'
params.bin_dir = '/proj/beyondfold/apps/alphafoldv2.3.1_pad'
params.af_dir = '/proj/beyondfold/apps/alphafold3'
params.af_flagfile = '/proj/beyondfold/apps/alphafold3/multimer.flag'
params.db_flagfile = '/proj/beyondfold/apps/alphafold3/databases.flag'
params.parse_flagfile = '/proj/beyondfold/apps/alphafold3/parse_features.flag'
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
    python ${params.bin_dir}/split_fasta.py $fasta split_fasta/
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
    clusterOptions "-A ${params.proj_id} -t 2:00:00 --gpus ${params.n_gpu} --reservation=devel"
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
        python ${params.bin_dir}/run_msa_tool.py $fasta mmseqs2 ${params.mmseqs_db} --out_dir ./ --gpu --mmseqs ${params.mmseqs_bin} --n_cpu $n_cpu $use_env
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
    python ${params.bin_dir}/af3/prepare_alignments.py $alignments AF_data/
    """
}

process parse_features {
    executor = 'local'
    publishDir {outputDir}, mode: 'copy'
    
    input:
    path fasta
    path af_data
    
    output:
    path "json_cache/*.json"
    
    script:
    """
    mkdir -p json_cache
    python ${params.bin_dir}/af3/parse_features.py --output_dir $af_data --fasta_paths $fasta --json_cache json_cache/ --flagfile ${params.parse_flagfile}
    """
}

process format_af_jobs {
    executor = 'local'
    publishDir {outputDir}

    input:
    path fasta
    path json_cache

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "AF_data_multimer/**.sh", emit: sh

    script:
    def file_list = params.file_list != '' ? "--file_list ${params.file_list}" : ''
    """
    python ${params.bin_dir}/af3/format_alphafold_jobs.py $fasta AF_data_multimer/ --conda_env ${params.conda_env} --json_dir ${outputDir}/json_cache --proj_id ${params.proj_id} --af3_path ${params.af_dir} $file_list --flagfiles ${params.af_flagfile} ${params.db_flagfile}
    """
}

process run_af_jobs {
    executor = "${params.executor}"
    clusterOptions "-A ${params.proj_id} -t 3-0 --gpus 1"
    publishDir {outputDir}
    
    input:
    path sbatch_script
    
    script:
    """
    if ${params.test}; then
        echo "Launching $sbatch_script"
    else
        sh $sbatch_script
    fi
    """
}

workflow {
    // align
    split_fasta_path = split_fasta(file(params.fasta))
    alignments_path = mmseqs_align(file(params.fasta))
    
    // convert
    af_data_path = convert_alignments(alignments_path)
    json_cache = parse_features(ln_fasta(split_fasta_path).flatten(), af_data_path).collect().flatten().take(1)
    
    // AF
    sbatch_scripts = format_af_jobs(split_fasta_path, json_cache).sh.collect().flatten()
    run_af_jobs(sbatch_scripts)
}
