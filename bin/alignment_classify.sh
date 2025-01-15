#!/bin/bash
SAMTOOLS_CONTAINER=/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/samtools/samtools%3A1.21--h50ea8bc_0

prefix=$1
bam=$2
seqid2taxid=$3
my_tax_ids=$4
nodes=$5
taxa_names=$6


# Determining child nodes based on parent tax ID provided
if [[ ! -z $nodes ]]; then

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

    tax_ids="$outfile"

fi

singularity exec $SAMTOOLS_CONTAINER samtools view -H $2 > "$prefix-classified.sam" #printing bam headers to output sam file

singularity exec $SAMTOOLS_CONTAINER samtools view $2 > "$prefix-temp.sam" #converting bam to sam for easier parsing in the loop below


declare -A taxa #creating dictionary to count the number of primary alignments present for each taxa 
num_reads=0

while IFS= read -r line; do #looping through each alignment in sam file

    flag=$(echo $line | awk '{print $2}') #grabbing value from "flag" column

    if [[ "$flag" == "0" ]]; then #confirm this is a primary alignment

        seq_id=$(echo $line | awk '{print $3}') #grabbing value of the reference the read aligned to
        tax_id=$(grep "$seq_id" $seqid2taxid | cut -f2 ) #converting reference seq id to tax id using seqid2taxid conversion file I stole from kraken2 

        if [[ ! -z $tax_id ]]; then

            ((num_reads++))

            if [[ $(grep -x "$tax_id" $tax_ids) ]]; then

                echo -e "$line" >> $prefix-classified.sam

            fi


            if [[ -v taxa["${tax_id}"] ]]; then

                ((taxa["$tax_id"]++))

            else
                echo $tax_id
                taxa[$tax_id]=1

            fi

        fi        

    else

        continue

    fi

done < "$prefix-temp.sam"

#converting sam file into bam file
singularity exec $SAMTOOLS_CONTAINER samtools view -@ 16 -S -b $prefix-classified.sam > $prefix-classified.bam


printing taxa count to a tsv file
>"$prefix-taxa-count.tsv"

for key in "${!taxa[@]}"; do
    percent=$(awk "BEGIN {print ${taxa["$key"]} / $num_reads}")
    echo -e "$key\t${taxa[$key]}\t$percent" >> $prefix-taxa-count.tsv
done

sort -k2,2nr $prefix-taxa-count.tsv > $prefix-taxa-count-by-ID-sorted.tsv
rm -f $prefix-taxa-count.tsv


#converting tax ids to taxa names
>$prefix-taxa-count-by-name-sorted.tsv

while IFS= read -r line; do

    tax_id=$(echo "$line" | awk '{print $1}')
    taxa_count=$(echo "$line" | awk '{print $2}')
    percent=$(echo "$line" | awk '{print $3}')
    taxa_name=$(grep "^$tax_id.[|].*scientific name" $taxa_names | cut -d "|" -f 2 | tr -d '\t')

    echo -e "$taxa_name\t$taxa_count\t$percent" >> $prefix-taxa-count-by-name-sorted.tsv   

done < $prefix-taxa-count-by-ID-sorted.tsv
