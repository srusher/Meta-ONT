process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/bbmap/bbmap%3A39.11--h92535d8_0' :
        'https://depot.galaxyproject.org/singularity/bbmap%3A39.06--h92535d8_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*deduped.fastq")      , optional:true, emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    reformat.sh in=$fastq out=${prefix}_filtered_again.fastq minlength=1

    dedupe.sh in=${prefix}_filtered_again.fastq out=${prefix}_filtered_deduped.fastq rmn=f
    """
}