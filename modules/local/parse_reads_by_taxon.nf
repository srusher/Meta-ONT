process PARSE_READS_BY_TAXON {
    tag "$meta.id"
    label 'process_medium'


    input:
    tuple val(meta), path(reads)
    path(tax_ids)

    output:
    tuple val(meta), path('*.fastq.gz') , optional:true, emit: taxon_reads

    script:
    def prefix = "${meta.id}"

    """
    
    bash "${projectDir}/scratch/parse_reads_by_taxon.sh" $prefix $reads $tax_ids

    """

}