include { convert_alignments ; parse_features_af3 ; collect_jsons ; format_jobs ; run_af3_jobs } from '../../../modules/local/af3'

workflow AF3 {
    take:
    alignments_path
    fasta_dir
    pair_list

    main:
    // convert
    af_data_path = convert_alignments(alignments_path)
    fastas = channel.fromPath("${fasta_dir}/*.fasta")
    jsons = parse_features_af3(fastas, af_data_path, params.af3_db_dir).collect()
    json_cache = collect_jsons(jsons)

    // AF
    jobs = format_jobs(fasta_dir, json_cache, pair_list, params.af3_db_dir).collect().flatten()
    run_af3_jobs(jobs, json_cache, params.af3_db_dir)
}
