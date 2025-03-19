process ALIGNMENT_CLASSIFY {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam), path(nanostats)
    path(seq2tax_map)
    val(filter_alignment_by_id)
    path(my_tax_ids)
    val(include_children)


    output:
    tuple val(meta), path('*-classified.bam') , optional:true, emit: classified_bam
    tuple val(meta), path('*-classified-plus-filtered.bam') , optional:true, emit: classified_plus_filtered_bam
    tuple val(meta), path('*ID-sorted.tsv') , optional:true, emit: id_tsv
    tuple val(meta), path('*classification-summary.tsv') , optional:true, emit: summary_tsv
    tuple val(meta), path('*alignment-classifiedreads.txt') , optional:true, emit: reads_tsv

    script:
    def prefix = "${meta.id}"

    """
    total_reads=\$(grep "Number of reads:" $nanostats | rev | awk '{print \$1}' | rev | cut -d '.' -f1 | tr -d ',')

    bash "${projectDir}/bin/alignment_classify.sh" $prefix $bam $seq2tax_map $filter_alignment_by_id $my_tax_ids $include_children ${params.ncbi_taxonomy_nodes} ${params.ncbi_taxonomy_names} \$total_reads ${projectDir} ${params.local_nodes_db}

    """

}