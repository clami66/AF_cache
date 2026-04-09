
process split_fasta {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta

    output:
    path "split_fasta/"

    script:
    """
    mkdir -p split_fasta/
    python ${params.af_cache_dir}/pipeline/common/split_fasta.py $fasta split_fasta/
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
    ln -s $fasta fasta
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
    
    output:
    path "alignments/"

    script:
    def use_env = params.use_env ? "--use-env" : ''
    def n_cpu = params.n_gpu > 0 ? params.n_cpus_per_gpu*params.n_gpu : params.max_cpus
    def use_gpu = params.n_gpu > 0 ? "--gpu " : ""
    """
    if ${params.test}; then
        mkdir -p alignments
        cp ${params.af_cache_dir}/test_data/alignments/*.a3m alignments/
    else
        python ${params.af_cache_dir}/pipeline/common/run_msa_tool.py $fasta \\
                                                                       mmseqs2 ${mmseqs_db}/ \\
                                                                       --out_dir ./ \\
                                                                       --mmseqs ${params.mmseqs_bin} \\
                                                                       --n_cpu $n_cpu \\
                                                                       $use_gpu \\
                                                                       $use_env
    fi
    """
}

