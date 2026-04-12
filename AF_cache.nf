#!/usr/bin/env nextflow

include { split_fasta ; ln_fasta ; mmseqs_align ; collect_pickles ; collect_jsons ; setup_mmseqs2_dbs } from './pipeline/common/modules'

process convert_alignments_af2 {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path alignments

    output:
    path "AF_data/"

    script:
    """
    python ${params.af_cache_dir}/pipeline/af2/prepare_alignments.py ${alignments} AF_data/
    """
}

process convert_alignments_af3 {
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
    python ${params.af_cache_dir}/pipeline/af2/parse_features.py --flagfile ${params.af2_flagfile} \\
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

process format_jobs_af2 {
    publishDir "${params.output_dir}"

    input:
    path fasta
    path pickle_cache
    path pair_list
    path af2_data_dir
    path template_mmcif_dir, stageAs: 'mmcif/*'
    path obsolete_pdbs_path, stageAs: 'obsolete/*'
    path pdb_seqres_database_path, stageAs: 'seqres/*'

    output:
    path "AF_data_multimer/", emit: 'dir'
    path "sbatch_scripts/**.sh", emit: sh

    script:
    def plist = pair_list.name != '.NO_FILE' ? "--file_list ${pair_list}" : ''
    def skip_templates = params.skip_templates ? "--notemplates" : ''
    """
    python ${params.af_cache_dir}/pipeline/af2/format_alphafold_jobs.py ${fasta} AF_data_multimer/ \\
                                                                        --pickle_dir ${pickle_cache} \\
                                                                        --write_fastas \\
                                                                        --af_path ${params.af_cache_dir} \\
                                                                        --flagfile ${params.af2_flagfile} \\
                                                                        --mmseqs2_uniref_database_path ${params.mmseqs_db}/uniref30_2302_db \\
                                                                        --mmseqs2_env_database_path ${params.mmseqs_db}/colabfold_envdb_202108_db \\
                                                                        --mmseqs2_binary_path ${params.mmseqs_bin} \\
                                                                        --template_mmcif_dir ${template_mmcif_dir} \\
                                                                        --obsolete_pdbs_path ${obsolete_pdbs_path} \\
                                                                        --pdb_seqres_database_path ${pdb_seqres_database_path} \\
                                                                        --data_dir ${af2_data_dir} \\
                                                                        ${skip_templates} \\
                                                                        ${plist}
    """
}

process format_jobs_af3 {
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
                                                                --af3_path ${params.af_cache_dir}/pipeline/af3 \\
                                                                --flagfile ${params.af3_flagfile} \\
                                                                --model_dir ${params.af3_model_dir} \\
                                                                ${plist}
    """
}

process run_af2_jobs {
    input:
    path sbatch_script
    path cache
    path template_mmcif_dir, stageAs: 'mmcif/*'
    path obsolete_pdbs_path, stageAs: 'obsolete/*'
    path pdb_seqres_database_path, stageAs: 'seqres/*'

    script:
    """
    sh ${sbatch_script}
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

process get_af2_params {
    output:
    path ".af2_params_ready"

    script:
    """
    if ! ls ${params.af2_data_dir}/params/params*.npz 1> /dev/null 2>&1; then
        ${params.af_cache_dir}/pipeline/af2/download_alphafold_params.sh ${params.af2_data_dir}
    fi

    touch .af2_params_ready
    """
}

workflow af2 {
    take:
    alignments_path
    fasta_links
    split_fasta_path
    pair_list

    main:
    // convert
    af_data_path = convert_alignments_af2(alignments_path)
    pickles = parse_features_af2(fasta_links, af_data_path, params.mmseqs_db, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).pkl.collect()
    pickle_cache = collect_pickles(pickles)

    // AF
    af2_data_dir = get_af2_params()
    sbatch_scripts = format_jobs_af2(split_fasta_path, pickle_cache, pair_list, af2_data_dir, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).sh.collect().flatten()
    run_af2_jobs(sbatch_scripts, pickle_cache, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path)
}

workflow af3 {
    take:
    alignments_path
    fasta_links
    split_fasta_path
    pair_list

    main:
    // convert
    af_data_path = convert_alignments_af3(alignments_path)
    jsons = parse_features_af3(fasta_links, af_data_path).collect()
    json_cache = collect_jsons(jsons)

    // AF
    sbatch_scripts = format_jobs_af3(split_fasta_path, json_cache, pair_list).sh.collect().flatten()
    run_af3_jobs(sbatch_scripts, json_cache)
}

workflow {
    // process inputs
    fasta = file(params.fasta)
    pair_list = file(params.pair_list)
    split_fasta_path = split_fasta(fasta)
    fasta_links = ln_fasta(split_fasta_path).flatten()    

    // align
    db_ready = setup_mmseqs2_dbs()
    alignments_path = mmseqs_align(fasta, params.mmseqs_db, db_ready)
    
    if ( params.af3 ) {
        af3(alignments_path, fasta_links, split_fasta_path, pair_list)
    } else {
        af2(alignments_path, fasta_links, split_fasta_path, pair_list)
    }
}
