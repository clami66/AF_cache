include { convert_alignments ; parse_features_af3 ; collect_jsons ; format_jobs ; run_af3_jobs } from '../../../modules/local/af3'
include { run_ipsae } from '../../../modules/local/ipsae'

workflow AF3 {
    take:
    alignments_path
    fasta_dir
    pair_list

    main:
    // convert
    af_data_path = convert_alignments(alignments_path)
    fastas = channel.fromPath("${fasta_dir}/*.fasta")
    jsons = parse_features_af3(fastas, af_data_path).collect()
    json_cache = collect_jsons(jsons)

    // AF
    jobs = format_jobs(fasta_dir, json_cache, pair_list).collect().flatten()
    predictions = run_af3_jobs(jobs, json_cache)

    // optional ipSAE scoring of the predicted complexes
    if (params.run_ipsae) {
        pair_dirs = predictions[0].flatten().collect()
        run_ipsae(pair_dirs, 'af3')
    }
}
