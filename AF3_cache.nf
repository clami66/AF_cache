#!/usr/bin/env nextflow

include { split_fasta ; ln_fasta ; mmseqs_align } from './pipeline/common/modules'

process convert_alignments {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path alignments

    output:
    path "AF_data/"

    script:
    """
    python ${params.af_cache_dir}/pipeline/af3/prepare_alignments.py ${alignments} AF_data/
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
    python ${params.af_cache_dir}/pipeline/af3/parse_features.py --output_dir ${af_data} \\
                                                        --fasta_paths ${fasta} \\
                                                        --json_cache json_cache/ \\
                                                        --flagfile ${params.af3_flagfile} \\
                                                        --undefok=num_diffusion_samples,model_dir
    """
}

process format_af_jobs {
    publishDir "${params.output_dir}"

    input:
    path fasta
    path json_cache
    path pair_list

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "sbatch_scripts/**.sh", emit: sh

    script:
    def plist = pair_list.name != '.NO_FILE' ? "--file_list ${pair_list}" : ''
    """
    python ${params.af_cache_dir}/pipeline/af3/format_alphafold_jobs.py ${fasta} AF_data_multimer/ \\
                                                                --json_dir ${json_cache} \\
                                                                --af3_path ${params.af3_install_path} \\
                                                                --flagfile ${params.af3_flagfile} \\
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

workflow {
    // align
    split_fasta_path = split_fasta(file(params.fasta))
    alignments_path = mmseqs_align(file(params.fasta), params.mmseqs_db)
    fasta_links = ln_fasta(split_fasta_path).flatten()
    // convert
    af_data_path = convert_alignments(alignments_path)
    jsons = parse_features_af3(fasta_links, af_data_path).collect()
    json_cache = collect_jsons(jsons)

    // AF
    sbatch_scripts = format_af_jobs(split_fasta_path, json_cache, file(params.pair_list)).sh.collect().flatten()
    run_af3_jobs(sbatch_scripts, json_cache)
}
