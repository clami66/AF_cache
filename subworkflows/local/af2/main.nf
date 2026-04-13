include { convert_alignments ; parse_features_af2 ; collect_pickles ; get_params ; format_jobs ; run_af2_jobs } from '../../../modules/local/af2'

workflow AF2 {
    take:
    alignments_path
    fasta_links
    split_fasta_path
    pair_list

    main:
    // convert
    af_data_path = convert_alignments(alignments_path)
    pickles = parse_features_af2(fasta_links, af_data_path, params.mmseqs_db, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).pkl.collect()
    pickle_cache = collect_pickles(pickles)

    // AF
    af2_data_ready = get_params()
    sbatch_scripts = format_jobs(split_fasta_path, pickle_cache, pair_list, af2_data_ready, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).sh.collect().flatten()
    run_af2_jobs(sbatch_scripts, pickle_cache, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path)
}
