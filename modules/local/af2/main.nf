process convert_alignments {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path alignments

    output:
    path "AF_data/"

    script:
    """
    af2_prepare_alignments.py ${alignments} AF_data/
    """
}

process parse_features_af2 {
    publishDir "${params.output_dir}", mode: 'copy'
    cache 'lenient'

    input:
    path fasta
    path af_data
    path mmseqs_db
    path template_mmcif_dir, stageAs: 'mmcif/*'
    path obsolete_pdbs_path, stageAs: 'obsolete/*'
    path pdb_seqres_database_path, stageAs: 'seqres/*'

    output:
    path "pickle_cache/**.pkl.gz", emit: pkl

    script:
    def skip_templates = params.skip_templates ? "--notemplates" : ''
    """
    # make sure that the custom installation of AF2.3 is found first
    mkdir -p pickle_cache
    af2_parse_features.py --flagfile ${params.af2_flagfile} \\
                        --output_dir ${af_data} \\
                        --fasta_paths ${fasta} \\
                        --mmseqs2_uniref_database_path ${mmseqs_db}/uniref30_2302_db \\
                        --mmseqs2_env_database_path ${mmseqs_db}/colabfold_envdb_202108_db \\
                        --template_mmcif_dir ${template_mmcif_dir} \\
                        --obsolete_pdbs_path ${obsolete_pdbs_path} \\
                        --pdb_seqres_database_path ${pdb_seqres_database_path} \\
                        ${skip_templates} \\
                        --undefok=data_dir,use_gpu_relax,models_to_relax,models_to_use,num_multimer_predictions_per_model,max_recycles \\
                        --pickle_cache pickle_cache/
    """
}

process format_jobs {

    input:
    path fasta
    path pickle_cache
    path pair_list
    path af2_data_ready
    path template_mmcif_dir, stageAs: 'mmcif/*'
    path obsolete_pdbs_path, stageAs: 'obsolete/*'
    path pdb_seqres_database_path, stageAs: 'seqres/*'

    output:
    path "chunk_*"

    script:
    def plist = pair_list ? "--file_list ${pair_list}" : ''
    def skip_templates = params.skip_templates ? "--notemplates" : ''
    def include_homomers = params.include_homomers ? '--include_homomers': ''
    """
    af2_format_jobs.py ${fasta} AF_data_multimer/ \\
                        --pickle_dir ${pickle_cache} \\
                        --write_fastas \\
                        --flagfile ${params.af2_flagfile} \\
                        --mmseqs2_uniref_database_path ${params.mmseqs_db}/uniref30_2302_db \\
                        --mmseqs2_env_database_path ${params.mmseqs_db}/colabfold_envdb_202108_db \\
                        --mmseqs2_binary_path ${params.mmseqs_bin} \\
                        --template_mmcif_dir ${template_mmcif_dir} \\
                        --obsolete_pdbs_path ${obsolete_pdbs_path} \\
                        --pdb_seqres_database_path ${pdb_seqres_database_path} \\
                        --data_dir ${params.af2_data_dir} \\
                        ${include_homomers} \\
                        ${skip_templates} \\
                        ${plist}
    """
}

process run_af2_jobs {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path chunk
    path cache
    path template_mmcif_dir, stageAs: 'mmcif/*'
    path obsolete_pdbs_path, stageAs: 'obsolete/*'
    path pdb_seqres_database_path, stageAs: 'seqres/*'

    output:
    path "AF_data_multimer/*"
    path "logs_AF2/${chunk}.log"

    script:
    """
    mkdir AF_data_multimer
    mkdir logs_AF2
    run_alphafold2.py --flagfile ${chunk}/chunk.flags --output_dir AF_data_multimer/ \$(cat ${chunk}/other.flags)
    cp .command.log logs_AF2/${chunk}.log
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

process get_params {
    output:
    path ".af2_params_ready"

    script:
    """
    if ! ls ${params.af2_data_dir}/params/params*.npz 1> /dev/null 2>&1; then
        download_alphafold_params.sh ${params.af2_data_dir}
    fi

    touch .af2_params_ready
    """
}
