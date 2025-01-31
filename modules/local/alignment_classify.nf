process ALIGNMENT_CLASSIFY {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam)
    path(seq2tax_map)
    path(my_tax_ids)
    val(include_children)


    output:
    tuple val(meta), path('*classified.bam') , optional:true, emit: bam
    tuple val(meta), path('*ID-sorted.tsv') , optional:true, emit: id_tsv
    tuple val(meta), path('*name-sorted.tsv') , optional:true, emit: name_tsv

    script:
    def prefix = "${meta.id}"

    """

    if [[ $include_children ]]; then

        bash "${projectDir}/bin/alignment_classify.sh" $prefix $bam $seq2tax_map $my_tax_ids ${params.ncbi_taxonomy_nodes} ${params.ncbi_taxonomy_names} ${params.strain_2_species}

    else

        bash "${projectDir}/bin/alignment_classify.sh" $prefix $bam $seq2tax_map $my_tax_ids "" ${params.ncbi_taxonomy_names}

    fi

    """

}