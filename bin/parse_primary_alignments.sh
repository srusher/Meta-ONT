#!/bin/bash

working_dir=$1
input_sam=$2
input_tax=$3
seqid2taxid=$4
filter_alignment_by_id=$5
tax_ids_i_want=$6
part=$7
prefix=$8
non_standard_reference=$9

declare -A taxa #creating dictionary to count the number of primary alignments present for each taxa 
num_classified_reads=0

>"$working_dir/$prefix-$part-classified.sam"
>"$working_dir/$prefix-$part-alignment-classifiedreads.txt"

if [[ "$filter_alignment_by_id" == "true" ]]; then

    >"$working_dir/$prefix-$part-classified-plus-filtered.sam"

fi

while IFS= read -r line; do #looping through each alignment in sam file

    flag=$(echo $line | awk '{print $2}') #grabbing value from "flag" column

    if [[ "$flag" == "0" || "$flag" == "16" ]]; then #confirm this is a primary alignment

        read_id=$(echo $line | awk '{print $1}')
        seq_id=$(echo $line | awk '{print $3}') #grabbing value of the reference the read aligned to

        if [[ $non_standard_reference == "true" ]]; then #if we're using a non-standard reference then there probably won't be a known sequence ID annotated in the fasta

            tax_id=$seq_id

        else

            tax_id=$(grep "$seq_id" $seqid2taxid | cut -f2 ) #converting reference seq id to tax id using a modified seqid2taxid conversion file I stole from kraken2 - modified by replacing all strain tax IDs with parent species tax IDs
        
        fi

        if [[ -n $tax_id ]]; then

            ((num_classified_reads++))

            echo "$line" >> $prefix-$part-classified.sam # DO NOT use '-e' flag here with echo! - This can cause CIGAR strings that contain the sequence "\n" to be split up into separate lines
            echo "$read_id\|$tax_id" >> $working_dir/$prefix-$part-alignment-classifiedreads.txt #adding read name and associated tax ID to separate file

            if [[ "$filter_alignment_by_id" == "true" ]]; then

                if [[ $(grep -x "$tax_id" $tax_ids_i_want) ]]; then # see if tax ID from this alignment is one of our specified tax IDs

                    echo "$line" >> $working_dir/$prefix-$part-classified-plus-filtered.sam

                fi

            fi                

            if [[ -v taxa["${tax_id}"] ]]; then

                ((taxa["$tax_id"]++))

            else
                
                taxa["$tax_id"]=1

            fi

        fi        

    else

        continue

    fi

done < "$input_sam"


# Print key-value pairs to the file
>"$working_dir/$prefix-$part-taxa-dictionary.tsv"

for key in "${!taxa[@]}"; do
    echo "$key|${taxa[$key]}" >> $working_dir/$prefix-$part-taxa-dictionary.tsv
done

# print classified reads count to file
echo $num_classified_reads > $working_dir/$prefix-$part-class-read-count.txt