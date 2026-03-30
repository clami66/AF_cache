#!/usr/bin/env nextflow

include { split_fasta; ln_fasta; mmseqs_align; } from './pipeline/common/modules'

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
    path mmseqs_db
    path template_mmcif_dir
    path obsolete_pdbs_path
    path pdb_seqres_database_path
    
    output:
    path "pickle_cache/**.pkl.gz", emit: pkl
    
    script:
    """
    # make sure that the custom installation of AF2.3 is found first
    mkdir -p pickle_cache
    python ${params.af_cache_dir}/pipeline/af2/parse_features.py --flagfile ${params.af2_flagfile} \\
                                                                 --output_dir $af_data \\
                                                                 --fasta_paths $fasta \\
                                                                 --mmseqs2_uniref_database_path ${mmseqs_db}/uniref30_2302_db \\
                                                                 --mmseqs2_env_database_path ${mmseqs_db}/colabfold_envdb_202108_db \\
                                                                 --template_mmcif_dir ${template_mmcif_dir} \\
                                                                 --obsolete_pdbs_path ${obsolete_pdbs_path} \\
                                                                 --pdb_seqres_database_path ${pdb_seqres_database_path} \\
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
    path template_mmcif_dir
    path obsolete_pdbs_path
    path pdb_seqres_database_path

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "AF_data_multimer/**.sh", emit: sh

    script:
    def file_list = params.file_list != '' ? "--file_list ${params.file_list}" : ''
    """
    python ${params.af_cache_dir}/pipeline/af2/format_alphafold_jobs.py $fasta AF_data_multimer/ \\
                                                                        --pickle_dir $pickle_cache \\
                                                                        --write_fastas \\
                                                                        --af_path ${params.af_cache_dir} \\
                                                                        --flagfile ${params.af2_flagfile} \\
                                                                        --mmseqs2_uniref_database_path ${params.mmseqs_db}/uniref30_2302_db \\
                                                                        --mmseqs2_env_database_path ${params.mmseqs_db}/colabfold_envdb_202108_db \\
                                                                        --mmseqs2_binary_path ${params.mmseqs_bin} \\
                                                                        --template_mmcif_dir ${template_mmcif_dir} \\
                                                                        --obsolete_pdbs_path ${obsolete_pdbs_path} \\
                                                                        --pdb_seqres_database_path ${pdb_seqres_database_path} \\
                                                                        --data_dir ${params.af2_data_dir} \\
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

process run_af2_jobs {
    executor = "${params.af_executor}"
    clusterOptions "${params.af_executor_flags}"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path sbatch_script
    path cache
    
    script:
    """
    sh $sbatch_script
    """
}

workflow {
    // align
    split_fasta_path = split_fasta(file(params.fasta))
    alignments_path = mmseqs_align(file(params.fasta), params.mmseqs_db)
    
    // convert
    af_data_path = convert_alignments(alignments_path)
    pickles = parse_features(ln_fasta(split_fasta_path).flatten(), af_data_path, params.mmseqs_db, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).pkl.collect()
    pickle_cache = collect_pickles(pickles)
    
    // AF
    sbatch_scripts = format_af_jobs(split_fasta_path, pickle_cache, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).sh.collect().flatten()
    run_af2_jobs(sbatch_scripts, pickle_cache)
}
