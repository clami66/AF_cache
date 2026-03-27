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
    python ${params.af_cache_dir}/pipeline/af2/prepare_alignments.py $alignments AF_data/
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
    path "pickle_cache/**.pkl.gz", emit: pkl
    
    script:
    """
    # make sure that the custom installation of AF2.3 is found first
    export PYTHONPATH="${params.af_cache_dir}:$PYTHONPATH"
    mkdir -p pickle_cache
    python ${params.af_cache_dir}/pipeline/af2/parse_features.py --flagfile ${params.af2_flagfile} \\
                                                                 --output_dir $af_data \\
                                                                 --fasta_paths $fasta \\
                                                                 --mmseqs2_uniref_database_path ${params.mmseqs_db}/uniref30_2302_db \\
                                                                 --mmseqs2_env_database_path ${params.mmseqs_db}/colabfold_envdb_202108_db \\
                                                                 --undefok=data_dir,use_gpu_relax,models_to_relax,models_to_use,num_multimer_predictions_per_model,max_recycles \\
                                                                 --pickle_cache pickle_cache/
    """
}

process format_af_jobs {
    executor = "${params.other_executor}"
    clusterOptions = "${params.other_executor_flags}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta
    path pickle_cache

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "AF_data_multimer/**.sh", emit: sh

    script:
    def file_list = params.file_list != '' ? "--file_list ${params.file_list}" : ''
    """
    python ${params.af_cache_dir}/pipeline/af2/format_alphafold_jobs.py $fasta AF_data_multimer/ \\
                                                                        --conda_env ${params.conda_env} \\
                                                                        --pickle_dir $pickle_cache \\
                                                                        --write_fastas \\
                                                                        --af_path ${params.af_cache_dir} \\
                                                                        --mmseqs2_uniref_database_path ${params.mmseqs_db}/uniref30_2302_db \\
                                                                        --mmseqs2_env_database_path ${params.mmseqs_db}/colabfold_envdb_202108_db \\
                                                                        --flagfile ${params.af2_flagfile} \\
                                                                        $file_list
    """
}

process collect_pickles {
    input:
    path "pickle_cache/*"

    output:
    path "pickle_cache"

    script:
    """
    echo "directory ready"
    """
}

workflow {
    // align
    split_fasta_path = split_fasta(file(params.fasta))
    alignments_path = mmseqs_align(file(params.fasta))
    
    // convert
    af_data_path = convert_alignments(alignments_path)
    pickles = parse_features(ln_fasta(split_fasta_path).flatten(), af_data_path).pkl.collect()
    pickle_cache = collect_pickles(pickles)
    
    // AF
    sbatch_scripts = format_af_jobs(split_fasta_path, pickle_cache).sh.collect().flatten()
    run_af_jobs(sbatch_scripts, pickle_cache)
}
