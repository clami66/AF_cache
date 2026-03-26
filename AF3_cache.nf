#!/usr/bin/env nextflow

outputDir = 'outputs'

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
    python ${params.af_cache_dir}/split_fasta.py $fasta split_fasta/
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
        python ${params.af_cache_dir}/run_msa_tool.py $fasta mmseqs2 ${params.mmseqs_db} --out_dir ./ --gpu --mmseqs ${params.mmseqs_bin} --n_cpu $n_cpu $use_env
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
    python ${params.af_cache_dir}/af3/prepare_alignments.py $alignments AF_data/
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
    conda run -p ${params.af3_conda_env} python ${params.af_cache_dir}/af3/parse_features.py --output_dir $af_data --fasta_paths $fasta --json_cache json_cache/ --flagfile ${params.af3_parse_flagfile}
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
    python ${params.af_cache_dir}/af3/format_alphafold_jobs.py $fasta AF_data_multimer/ --json_dir ${outputDir}/json_cache --proj_id ${params.proj_id} --af3_path ${params.af3_dir} $file_list --flagfiles ${params.af3_flagfile} ${params.af3_db_flagfile}
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
