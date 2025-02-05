process KRAKEN_ALIGNMENT_COMPARISON {
    tag "$meta.id"
    label 'process_medium'
    //errorStrategy 'ignore'

    input:
    tuple val(meta), path(kraken_report), path(alignment_report)

    output:
    tuple val(meta), path('*.png') , optional:true, emit: plot

    script:
    def prefix = "${meta.id}"

    """

    bash "${projectDir}/bin/alignment_kraken_comparison.sh" $prefix $kraken_report $alignment_report "${projectDir}/bin/kraken_alignment_plot.R"

    """

}