#!/bin/bash
SAMTOOLS_CONTAINER="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/samtools/samtools%3A1.21--h50ea8bc_0"
SQLITE3_CONTAINER="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/sqlite3/sqlite%3A3"

prefix=$1
bam=$2
seqid2taxid=$3
filter_alignment_by_id=$4
my_tax_ids=$5
include_children=$6
nodes=$7
taxa_names=$8
total_reads=$9
project_dir=${10}
nodes_sqlite=${11}

# Determining child nodes based on parent tax ID provided
if [[ "$filter_alignment_by_id" == "true" && "$include_children" == "true" ]]; then

    echo "Collecting all children tax ids from SQL database"

    >"$prefix-parent_and_child_ids.txt"

    for i in $(cat $my_tax_ids); do

        data=$(singularity exec --bind /scicomp $SQLITE3_CONTAINER sqlite3 $nodes_sqlite "SELECT * FROM TAX_IDS WHERE parent_id = "$i"")

        parent_id="$i"        
        echo $parent_id >> $prefix-parent_and_child_ids.txt

        if [[ -n $data ]]; then

            child_ids=$(echo $data | cut -d '|' -f3 | sed 's/,/ /g' )

            for i in $child_ids; do

                echo $i >> $prefix-parent_and_child_ids.txt
            
            done
        
        fi

    done

    tax_ids_i_want="$prefix-parent_and_child_ids.txt"

else

    tax_ids_i_want=""

fi

echo "Parsing BAM headers"

singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view -H $2 > "$prefix-classified.sam" #printing bam headers to output sam file

echo "Converting BAM to SAM with no headers"

singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view $2 > "$prefix-temp.sam" #converting bam to sam for easier parsing in the loop below

if [[ "$filter_alignment_by_id" == "true" ]]; then

    singularity exec --bind /scicomp $SAMTOOLS_CONTAINER samtools view -H $2 > "$prefix-classified-plus-filtered.sam" #printing bam headers to output sam file

fi

chunks=10

# Get the total number of lines in the sam file
total_lines=$(wc -l $prefix-temp.sam | cut -d ' ' -f 1)

# Calculate the number of lines per chunk (divided into n chunks)
div=$((total_lines/chunks))
start=1
fin=$(($div+1))

function parse_sam() {

    for i in $(seq 1 $chunks); do

        working_dir=$(pwd)
        sam_chunk="$working_dir/sam_chunk_p$i"
        tax_count="$working_dir/tax_count_p$i"
        tax_ids_i_want="$working_dir/$tax_ids_i_want"

        if [ "$i" -eq $chunks ]; then
        
            #sometimes the SAM file cannot be broken up into n even chunks - so on the nth chunk we need to ensure we include all of the remaining lines in the file
            sed -n ''"$start"','"$total_lines"'p' $prefix-temp.sam > sam_chunk_p$i

        else 
            # Extract lines from the input file into the output file
            sed -n ''"$start"','"$fin"'p' $prefix-temp.sam > sam_chunk_p$i

        fi

        session_name="parse_sam_$i-$prefix"

        tmux new-session -d -s $session_name "bash $project_dir/bin/parse_primary_alignments.sh $working_dir $sam_chunk $tax_count $seqid2taxid $filter_alignment_by_id $tax_ids_i_want part-$i $prefix"         
    
        # Update the start and finish line numbers for the next chunk
        start=$((start + div + 1))
        fin=$((fin + div))

    done
}

echo "Splitting up SAM parsing amongst tmux sessions"

#execute above function
parse_sam

# only continue script once all tmux sessions have finished processing their respective chunks
continue_code="false"
while [[ $continue_code == "false" ]]; do

    if [[ $(tmux list-sessions | cut -d ':' -f1 | grep "parse_sam" | wc -l) -gt 0 ]]; then

        sleep 5
    
    else

        continue_code="true"
    
    fi

done

echo "Concatenating respective tmux output files"

#combining tmux output files into one file
cat ./*"-part"*"alignment-classifiedreads.txt" > $prefix-alignment-classifiedreads.txt
cat ./*"-part"*"classified.sam" >> $prefix-classified.sam

if [[ "$filter_alignment_by_id" == "true" ]]; then

    cat ./*"-part"*"classified-plus-filtered.sam" >> $prefix-classified-plus-filtered.sam

fi


# calculating total classified reads from tmux output files
cat ./*"-part"*"class-read-count.txt" > $prefix-class-read-count.txt
num_classified_reads=0
while IFS= read -r line; do

    num_classified_reads=$((num_classified_reads + line))

done < "$prefix-class-read-count.txt"

# creating aggregate taxa dict from tmux output files
cat ./*"-part"*"taxa-dictionary.tsv" > $prefix-taxa-dictionary-all.tsv
declare -A taxa

while IFS= read -r line; do

    tax_id=$(echo $line | cut -d '|' -f1)
    count=$(echo $line | cut -d '|' -f2)

    if [[ -v taxa["${tax_id}"] ]]; then

        value=${taxa[$tax_id]}

        echo "value: $value"

        total=$((value + count))
        taxa["$tax_id"]="$total"

    else
        
        taxa["$tax_id"]=$count

    fi

done < "$prefix-taxa-dictionary-all.tsv"



if [[ "$filter_alignment_by_id" == "true" ]]; then

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

# cleaning up temporary files
rm -f $prefix-part-* 
rm -f sam_chunk*
rm -f $prefix-*.sam