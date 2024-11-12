process PARSE_READS_BY_TAXON {
    tag "$meta.id"
    label 'process_medium'


    input:
    tuple val(meta), path(report), path(reads)

    output:
    path('*.fastq') , emit: taxon_reads
    path('*.fasta') , emit: taxon_fasta
    path('*_revised.txt') , emit: blast_results

    script:
    def prefix = "${meta.id}"
    """
    cat ${report} | awk '\$2 == "5755" || \$2 == "1257118" || \$2 == "9606"' | awk '{print \$1}' > $prefix-filtered-taxon-ids

    >$prefix-taxon-filtered-reads.fastq

    for i in \$(cat $prefix-filtered-taxon-ids); do zcat ${reads} | awk "/@\$i/ {count=3; print; next} count > 0 {print; count--}" >> $prefix-taxon-filtered-reads.fastq; done

    cat $prefix-taxon-filtered-reads.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $prefix-taxon-filtered-reads.fasta

    if [ ! -s $prefix-taxon-filtered-reads.fasta ]; then

        echo "no taxon reads found"

    else

        BLASTDB=/scicomp/reference/ncbi-blast-taxdb
        blast_args="6 qseqid sacc pident length mismatch evalue bitscore stitle"
        singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/blast/blast-2.15.0--pl5321h6f7f691_1 blastn \\
            -db "${params.blast_db}" \\
            -query $prefix-taxon-filtered-reads.fasta \\
            -outfmt "\$blast_args" \\
            -num_threads 16 \\
            -evalue ${params.blast_evalue} -perc_identity ${params.blast_perc_identity} -max_target_seqs ${params.blast_target_seqs} \\
            -out ${prefix}.txt

        mod_blast_args="\${blast_args/6 /}"
        mod_blast_args=\$(echo "\$mod_blast_args" | sed 's/ /\t/g')
        echo -e "\$mod_blast_args" > ${prefix}_revised.txt

        cat ${prefix}.txt >> ${prefix}_revised.txt

        rm -f ${prefix}.txt
    
    fi
    """

}