process setup_mmseqs2_dbs {
    output:
    path ".databases_ready"

    script:
    def gpu_var = params.mmseqs_n_gpu > 0 ? 'GPU=1' : 'GPU=0'
    def template_var = params.skip_templates ? 'SKIP_TEMPLATES=1' : ''
    """
    if [[ ! -e "${params.mmseqs_db}" ]]; then
        mkdir -p "${params.mmseqs_db}"
    fi

    if [[ ! -e "${params.mmseqs_db}/PDB_MMCIF_READY" ]] || [[ ! -e "${params.mmseqs_db}/UNIREF30_READY" ]] || [[ ! -e "${params.mmseqs_db}/COLABDB_READY" ]]; then
        ${gpu_var} ${template_var} FAST_PREBUILT_DATABASES=1 ${params.af_cache_dir}/pipeline/common/setup_databases.sh ${params.mmseqs_db}
    fi

    touch .databases_ready
    """
}

process split_fasta {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta

    output:
    path "split_fasta/"

    script:
    """
    mkdir -p split_fasta/
    python ${params.af_cache_dir}/pipeline/common/split_fasta.py ${fasta} split_fasta/
    """
}

process ln_fasta {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta

    output:
    path "fasta/*.fasta"

    script:
    """
    ln -s ${fasta} fasta
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

process collect_jsons {
    input:
    path "json_cache/*"

    output:
    path "json_cache"

    script:
    """
    echo "directory ready"
    """
}

process mmseqs_align {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta
    path mmseqs_db
    path db_ready

    output:
    path "alignments/"

    script:
    def use_env = params.mmseqs_use_env ? "--use-env" : ''
    def n_cpu = params.mmseqs_max_cpus
    def use_gpu = params.mmseqs_n_gpu > 0 ? "--gpu " : ""
    """
    if ${params.test}; then
        mkdir -p alignments
        cp ${params.af_cache_dir}/test_data/alignments/*.a3m alignments/
    else
        python ${params.af_cache_dir}/pipeline/common/run_msa_tool.py ${fasta} \\
                                                                       mmseqs2 ${mmseqs_db}/ \\
                                                                       --out_dir ./ \\
                                                                       --n_cpu ${n_cpu} \\
                                                                       ${use_gpu} \\
                                                                       ${use_env}
    fi
    """
}
