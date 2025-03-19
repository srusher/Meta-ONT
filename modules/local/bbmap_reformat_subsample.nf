process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/bbmap/bbmap%3A39.11--h92535d8_0' :
        'https://depot.galaxyproject.org/singularity/bbmap%3A39.06--h92535d8_0' }"

    input:
    tuple val(meta), path(fastq)
    val(num_subsamples)

    output:
    tuple val(meta), path("*_subsample*.fastq.gz*")      , optional:true, emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    reformat.sh in=${fastq[0]} \
        out=${prefix}_subsample-${num_subsamples}.fastq \
        samplereadstarget=${num_subsamples} \
        sampleseed=13

    gzip -f *_subsample*.fastq
    """
}