#!/bin/bash

SAMTOOLS_CONTAINER=/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/samtools/samtools%3A1.21--h50ea8bc_0

prefix=$1
bam=$2
seqid2taxid=$3
tax_ids=$4

singularity exec $SAMTOOLS_CONTAINER samtools view -H $2 > "$prefix-classified.sam"

sam=$(singularity exec $SAMTOOLS_CONTAINER samtools view $2)

while IFS= read -r line; do

    seq_id=$(echo $line | awk '{print $3}')
    tax_id=$(grep "$seq_id" $seqid2taxid | cut -f2 )

    if [[ $(grep -x "$tax_id" $tax_ids) ]]; then

        echo -e "$line" >> $prefix-classified.sam

    fi

done <<< "$sam"

singularity exec $SAMTOOLS_CONTAINER samtools view -@ 16 -S -b $prefix-classified.sam > $prefix-classified.bam