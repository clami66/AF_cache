#!/usr/bin/env nextflow

include { split_fasta; ln_fasta; mmseqs_align; run_af_jobs } from './pipeline/common/modules'

process convert_alignments {
    executor = "${params.other_executor}"
    clusterOptions = "${params.other_executor_flags}"
    publishDir "${params.output_dir}", mode: 'copy'

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
    executor = "${params.other_executor}"
    clusterOptions = "${params.other_executor_flags}"
    publishDir "${params.output_dir}", mode: 'copy'
    
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
    executor = "${params.other_executor}"
    clusterOptions = "${params.other_executor_flags}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta
    path json_cache

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "AF_data_multimer/**.sh", emit: sh

    script:
    def file_list = params.file_list != '' ? "--file_list ${params.file_list}" : ''
    """
    python ${params.af_cache_dir}/af3/format_alphafold_jobs.py $fasta AF_data_multimer/ \\
                                                                --json_dir ${params.output_dir}/json_cache \\
                                                                --af3_path ${params.af3_dir} $file_list \\
                                                                --flagfiles ${params.af3_flagfile} \\
                                                                ${params.af3_db_flagfile}
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
