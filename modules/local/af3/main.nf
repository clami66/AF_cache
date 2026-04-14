process convert_alignments {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path alignments

    output:
    path "AF_data/"

    script:
    """
    af3_prepare_alignments.py ${alignments} AF_data/
    """
}

process parse_features_af3 {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta
    path af_data

    output:
    path "json_cache/*.json"

    script:
    """
    mkdir -p json_cache
    af3_parse_features.py --output_dir ${af_data} \\
                        --fasta_paths ${fasta} \\
                        --json_cache json_cache/ \\
                        --flagfile ${params.af3_flagfile} \\
                        --undefok=num_diffusion_samples,model_dir
    """
}

process format_jobs {
    publishDir "${params.output_dir}"

    input:
    path fasta
    path json_cache
    path pair_list

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "sbatch_scripts/**.sh", emit: sh

    script:
    def plist = pair_list ? "--file_list ${pair_list}" : ''
    """
    af3_format_jobs.py ${fasta} AF_data_multimer/ \\
                        --json_dir ${json_cache} \\
                        --flagfile ${params.af3_flagfile} \\
                        --model_dir ${params.af3_model_dir} \\
                        --include_homomers ${params.include_homomers} \\
                        ${plist}
    """
}

process collect_jsons {
    input:
    path "json_cache/*"

    output:
    path "json_cache"

    script:
    """
    echo "directory ready"
    """
}

process run_af3_jobs {
    input:
    path sbatch_script
    path cache

    script:
    """
    sh ${sbatch_script}
    """
}
