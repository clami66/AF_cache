include { convert_alignments ; parse_features_af2 ; collect_pickles ; get_params ; format_jobs ; run_af2_jobs } from '../../../modules/local/af2'
include { run_ipsae } from '../../../modules/local/ipsae'

workflow AF2 {
    take:
    alignments_path
    fasta_dir
    pair_list

    main:
    // convert
    af_data_path = convert_alignments(alignments_path)
    fastas = channel.fromPath("${fasta_dir}/*.fasta")
    pickles = parse_features_af2(fastas, af_data_path, params.mmseqs_db, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).pkl.collect()
    pickle_cache = collect_pickles(pickles)

    // AF
    af2_data_ready = get_params()
    jobs = format_jobs(fasta_dir, pickle_cache, pair_list, af2_data_ready, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path).collect().flatten()
    predictions = run_af2_jobs(jobs, pickle_cache, params.template_mmcif_dir, params.obsolete_pdbs_path, params.pdb_seqres_database_path)

    // optional ipSAE scoring of the predicted complexes
    if (params.run_ipsae) {
        pair_dirs = predictions[0].flatten().collect()
        run_ipsae(pair_dirs, 'af2')
    }
}
