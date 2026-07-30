process run_ipsae {
    publishDir "${params.output_dir}/metrics", mode: 'copy'

    input:
    path predictions, stageAs: 'predictions/*'
    val backend

    output:
    path "ipsae.tsv"

    script:
    """
    run_ipsae.py \\
        --predictions-dir predictions/ \\
        --backend ${backend} \\
        --ipsae-script ${params.ipsae_script} \\
        --pae-cutoff ${params.ipsae_pae_cutoff} \\
        --dist-cutoff ${params.ipsae_dist_cutoff} \\
        --workers ${task.cpus} \\
        --output ipsae.tsv
    """
}
