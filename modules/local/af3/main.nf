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
    path af3_db_dir

    output:
    path "json_cache/*.json"

    script:
    def skip_templates = params.skip_templates ? "--notemplates" : ''
    """
    mkdir -p json_cache
    af3_parse_features.py --output_dir ${af_data} \\
                        --fasta_paths ${fasta} \\
                        --json_cache json_cache/ \\
                        --flagfile ${params.af3_flagfile} \\
                        --db_dir ${af3_db_dir} \\
                        ${skip_templates} \\
                        --undefok=num_diffusion_samples,model_dir
    """
}

process format_jobs {

    input:
    path fasta
    path json_cache
    path pair_list
    path af3_db_dir

    output:
    path "chunk_*"

    script:
    def plist = pair_list ? "--file_list ${pair_list}" : ''
    def include_homomers = params.include_homomers ? '--include_homomers': ''
    """
    af3_format_jobs.py ${fasta} AF_data_multimer/ \\
                        --json_dir ${json_cache} \\
                        --flagfile ${params.af3_flagfile} \\
                        --model_dir ${params.af3_model_dir} \\
                        --n_seeds ${params.af3_n_prediction_seeds} \\
                        --db_dir ${af3_db_dir} \\
                        ${include_homomers} \\
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
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path chunk
    path cache
    path af3_db_dir

    output:
    path "AF3_data_multimer/*"
    path "logs_AF3/${chunk}.log"

    script:
    """
    mkdir AF3_data_multimer
    mkdir logs_AF3
    run_alphafold3.py --flagfile ${chunk}/chunk.flags --output_dir AF3_data_multimer/ \$(cat ${chunk}/other.flags)
    cp .command.log logs_AF3/${chunk}.log
    """
}
