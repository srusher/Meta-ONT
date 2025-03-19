process READS_TO_BLASTDB {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(fastq)
    path  gene

    output:
    tuple val(meta), path('*sorted.txt'), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gwa_parent="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL"

    awk 'NR%4==1 {print ">" substr(\$0, 2)} NR%4==2 {print}' $fastq > ${prefix}_kraken-translate_filtered.fasta

    singularity exec --bind /scicomp \$gwa_parent/singularity/blast/blast-2.15.0--pl5321h6f7f691_1 makeblastdb -in ${prefix}_kraken-translate_filtered.fasta -dbtype nucl -out nucl_reads

    BLAST_CONTAINER="/apps/standalone/singularity/blast/blast-2.15.0--pl5321h6f7f691_1"
    BLASTDB=/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/Projects/Long_Read_Analysis/data/blast/taxonomy
    blast_args="6 qseqid qstart qend sseqid sstart send pident length mismatch evalue bitscore sacc stitle staxids"

    singularity exec --bind /scicomp \$BLAST_CONTAINER blastn -db "./nucl_reads" -query $gene -outfmt "\$blast_args" -num_threads 16 -out ${prefix}.txt

    mod_blast_args="\${blast_args/6 /}"
    mod_blast_args=\$(echo "\$mod_blast_args" | sed 's/ /\t/g')
    echo -e "\$mod_blast_args" > ${prefix}_blast-hits_sorted.txt

    awk '{print \$0}' ${prefix}.txt | sort -k11,11 -nr >> ${prefix}_blast-hits_sorted.txt

    rm -f ${prefix}.txt
    rm -f ${prefix}_revised.txt
    """
}
