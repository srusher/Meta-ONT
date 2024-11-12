process RENAME_KRAKEN_READS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq*"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bash_reads=$reads
    mod_reads="\${bash_reads/.unclassified/_filtered}"
    mod_reads="\${mod_reads/.classified/_filtered}"
    mv $reads \$mod_reads
    """

}
