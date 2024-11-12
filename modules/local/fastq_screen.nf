process FASTQSCREEN {
    tag "$meta.id"
    label 'process_medium'


    input:
    tuple val(meta), path(reads)
    path(conf)

    output:
    path('*') , emit: fastq_screen_report

    script:
    def prefix = "${meta.id}"
    """
    singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/fastq-screen/custom/fastq_screen-0.15.3_modified.sif fastq_screen ${reads} --conf ./${conf} --aligner minimap2 --tag ${reads}
    """

}