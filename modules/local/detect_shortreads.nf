process DETECT_SHORTREADS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    val(meta.short_1), emit: detection
    path('*.txt') , emit: file1

    when:
    "${meta.short_1}" == "*fastq*"

    script:
    """
    echo "${meta.short_1}" > ${meta.id}_detection.txt
    """

}
