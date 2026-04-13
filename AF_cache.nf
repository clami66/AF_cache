#!/usr/bin/env nextflow

include { AF2 } from './subworkflows/local/af2'
include { AF3 } from './subworkflows/local/af3'
include { split_fasta ; ln_fasta ; mmseqs_align ; setup_mmseqs2_dbs } from './modules/local/common'

workflow {
    // process inputs
    fasta = file(params.fasta)
    pair_list = params.pair_list ? file(params.pair_list) : params.pair_list
    split_fasta_path = split_fasta(fasta)
    fasta_links = ln_fasta(split_fasta_path).flatten()

    // align
    db_ready = setup_mmseqs2_dbs()
    alignments_path = mmseqs_align(fasta, params.mmseqs_db, db_ready)

    if (params.af3) {
        AF3(alignments_path, fasta_links, split_fasta_path, pair_list)
    }
    else {
        AF2(alignments_path, fasta_links, split_fasta_path, pair_list)
    }
}
