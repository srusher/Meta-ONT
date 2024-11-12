process ALIGNMENT_CLASSIFY {
    tag "$meta.id"
    label 'process_high'
    //errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam)
    path(seq2tax_map)
    path(tax_ids)

    output:
    tuple val(meta), path('*classified.bam') , optional:true, emit: bam

    script:
    def prefix = "${meta.id}"

    """

    bash "${projectDir}/bin/alignment_classify.sh" $prefix $bam $seq2tax_map $tax_ids

    """

}