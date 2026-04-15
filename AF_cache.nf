#!/usr/bin/env nextflow

include { AF2 } from './subworkflows/local/af2'
include { AF3 } from './subworkflows/local/af3'
include { merge_fastas ; mmseqs_align ; setup_mmseqs2_dbs } from './modules/local/common'

workflow {
    // process inputs
    fasta_dir = file(params.fasta)
    merged_fasta = merge_fastas(fasta_dir)
    pair_list = params.pair_list ? file(params.pair_list) : params.pair_list

    // align
    db_ready = setup_mmseqs2_dbs()
    alignments_path = mmseqs_align(merged_fasta, params.mmseqs_db, db_ready)

    if (params.af3) {
        AF3(alignments_path, fasta_dir, pair_list)
    }
    else {
        AF2(alignments_path, fasta_dir, pair_list)
    }
}
