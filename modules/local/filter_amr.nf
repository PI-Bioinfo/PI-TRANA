process DIAMOND_FILTER {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(diamond_out)
    
    output:
    tuple val(meta), path("*_amr_filter.tsv"), emit: amr
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def basename = diamond_out.baseName  // gets filename without extension
    """
    awk '\$3>=70 && \$4>=30 && \$11<=1e-5 && \$12>=50' ${diamond_out} > ${basename}_amr_filter.tsv
    """
}