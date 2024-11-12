#!/bin/bash

report=$1
prefix=$2
reads=$3

cat $report | awk '$2 == "211522"' | awk '{print $1}' > $prefix-filtered-taxon-ids

>"$prefix-taxon-filtered-reads.fastq"

for i in $(cat $prefix-filtered-taxon-ids); do 
    
    zcat $reads | awk "/@$i/ {count=3; print; next} count > 0 {print; count--}" >> $prefix-taxon-filtered-reads.fastq; 
    
done

cat $prefix-taxon-filtered-reads.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $prefix-taxon-filtered-reads.fasta
