#!/bin/bash
SAMTOOLS_CONTAINER="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/samtools/samtools%3A1.21--h50ea8bc_0"

prefix=$1
bam=$2
seqid2taxid=$3
filter_alignment_by_id=$4
my_tax_ids=$5
include_children=$6
nodes=$7
taxa_names=$8
total_reads=$9

# Determining child nodes based on parent tax ID provided
if [[ "$filter_alignment_by_id" == "true" && "$include_children" == "true" ]]; then

    outfile="$prefix-parent_and_child_ids.txt"
    >$outfile

    while IFS= read -r id; do #iterating through tax ids listed in the tax id input file

        echo "$id" >> $outfile
        tax_array=("$id") #creating array for tax ids
        loop_again=true

        while $loop_again; do

            count=0

            for i in "${tax_array[@]}"; do

                if [[ $count -eq 0 ]]; then

                    tax_array=() #clearing out array to rebuild it with child tax ids for the next for loop iteration
                    ((count++))

                fi

                child_ids="$(awk -v id="$i" '$3 == id { print $1 }' $nodes)" #finding child nodes that have the "parent tax id" column set to current tax id

                if [[ ! -z $child_ids ]]; then #checking to see if we found any child nodes

                    for x in $child_ids; do

                        echo "$x" >> $outfile
                        tax_array+=("$x") #appending child tax ids to array to use in 

                    done                        

                else 

                    continue

                fi

            done

            if [ ${#tax_array[@]} -eq 0 ]; then

                loop_again=false

            fi

        done

    done < "$my_tax_ids"

    tax_ids_i_want="$outfile"

fi

singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view -H $2 > "$prefix-classified.sam" #printing bam headers to output sam file

singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view $2 > "$prefix-temp.sam" #converting bam to sam for easier parsing in the loop below

if [[ "$filter_alignment_by_id" == "true" ]]; then

    singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view -H $2 > "$prefix-classified-plus-filtered.sam" #printing bam headers to output sam file

fi

>"$prefix-alignment-classifiedreads.txt"

declare -A taxa #creating dictionary to count the number of primary alignments present for each taxa 
num_classified_reads=0

while IFS= read -r line; do #looping through each alignment in sam file

    flag=$(echo $line | awk '{print $2}') #grabbing value from "flag" column

    if [[ "$flag" == "0" || "$flag" == "16" ]]; then #confirm this is a primary alignment

        read_id=$(echo $line | awk '{print $1}')
        seq_id=$(echo $line | awk '{print $3}') #grabbing value of the reference the read aligned to
        tax_id=$(grep "$seq_id" $seqid2taxid | cut -f2 ) #converting reference seq id to tax id using a modified seqid2taxid conversion file I stole from kraken2 - modified by replacing all strain tax IDs with parent species tax IDs

        if [[ -n $tax_id ]]; then

            ((num_classified_reads++))

            echo "$line" >> $prefix-classified.sam # DO NOT use '-e' flag here with echo! - This can cause CIGAR strings that contain the sequence "\n" to be split up into separate lines
            echo -e "$read_id\t$tax_id" >> $prefix-alignment-classifiedreads.txt #adding read name and associated tax ID to separate file

            if [[ "$filter_alignment_by_id" == "true" ]]; then

                if [[ $(grep -x "$tax_id" $tax_ids_i_want) ]]; then # see if tax ID from this alignment is one of our specified tax IDs

                    echo "$line" >> $prefix-classified-plus-filtered.sam

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

done < "$prefix-temp.sam"


if [[ "$filter_alignment_by_id" == "true" ]]; then # if we want to filter the new sam file for certain tax IDs proceed with this block

    #converting sam file into bam file
    singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view -@ 16 -S -b $prefix-classified-plus-filtered.sam > $prefix-classified-plus-filtered.bam

fi

#converting sam file into bam file
singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view -@ 16 -S -b $prefix-classified.sam > $prefix-classified.bam


#printing taxa count to a tsv file
>"$prefix-taxa-count.tsv"

for key in "${!taxa[@]}"; do
    percent=$(awk "BEGIN {print ${taxa["$key"]} / $total_reads}")
    echo -e "$key\t${taxa[$key]}\t$percent" >> $prefix-taxa-count.tsv
done

sort -k2,2nr $prefix-taxa-count.tsv > $prefix-taxa-count-by-ID-sorted.tsv
rm -f $prefix-taxa-count.tsv


#converting tax ids to taxa names
>$prefix-alignment-classification-summary.tsv

unclassified_reads=$((total_reads - num_classified_reads))
unclassified_reads_percent=$(awk "BEGIN {print ($unclassified_reads / $total_reads) * 100}")
echo -e "unclassified\t$unclassified_reads\t$unclassified_reads_percent" >> $prefix-alignment-classification-summary.tsv

while IFS= read -r line; do

    tax_id=$(echo "$line" | awk '{print $1}')
    taxa_count=$(echo "$line" | awk '{print $2}')
    percent=$(echo "$line" | awk '{print $3}' | awk '{$1=$1*100; print $0}')
    taxa_name=$(grep "^$tax_id.[|].*scientific name" $taxa_names | cut -d "|" -f 2 | tr -d '\t') # this line is grabbing the scientific name from the names.dmp taxonomy file

    sed -i "s|$tax_id|$taxa_name|g" $prefix-alignment-classifiedreads.txt # replacing tax id with scientific name in the classified reads file
    echo -e "$taxa_name\t$taxa_count\t$percent" >> $prefix-alignment-classification-summary.tsv   

done < $prefix-taxa-count-by-ID-sorted.tsv