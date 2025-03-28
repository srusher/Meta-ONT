process NANOPLOT {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*report.html")                , emit: html
    tuple val(meta), path("*dot.png") , optional: true, emit: dot
    tuple val(meta), path("*kde.png") , optional: true, emit: kde
    tuple val(meta), path("*.txt")                 , emit: txt
    //tuple val(meta), path("*.log")                 , emit: log
    //path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_file = ("$ontfile".endsWith(".fastq.gz")) ? "--fastq ${ontfile}" : ("$ontfile".endsWith(".fastq")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        --prefix "${prefix}_" \\
        $input_file

    if [[ "$input_file" == *"trimmed"* ]]; then

        mv ${prefix}_NanoStats.txt ${prefix}_NanoStats_trimmed.txt

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch LengthvsQualityScatterPlot_dot.html
    touch LengthvsQualityScatterPlot_kde.html
    touch NanoPlot-report.html
    touch NanoPlot_20240301_1130.log
    touch NanoStats.txt
    touch Non_weightedHistogramReadlength.html
    touch Non_weightedLogTransformed_HistogramReadlength.html
    touch WeightedHistogramReadlength.html
    touch WeightedLogTransformed_HistogramReadlength.html
    touch Yield_By_Length.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
