process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast%3A2.15.0--pl5321h6f7f691_0':
        'biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    BLASTDB=/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/Projects/Long_Read_Analysis/data/blast/taxonomy
    blast_args="6 qseqid qstart qend sseqid sstart send pident length mismatch evalue bitscore sacc stitle staxids"
    blastn \\
        -db "${params.blast_db}" \\
        -query $fasta \\
        -outfmt "\$blast_args" \\
        -num_threads 16 \\
        $args \\
        -out ${prefix}.txt

    mod_blast_args="\${blast_args/6 /}"
    mod_blast_args=\$(echo "\$mod_blast_args" | sed 's/ /\t/g')
    echo -e "\$mod_blast_args" > ${prefix}_blast-hits_sorted.txt

    awk '{print \$0}' ${prefix}.txt | sort -k11,11 -nr >> ${prefix}_blast-hits_sorted.txt

    rm -f ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
