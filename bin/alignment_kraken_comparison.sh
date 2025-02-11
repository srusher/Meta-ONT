#!/bin/bash
R_CONTAINER="/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/R/ggplot_dplyr/ggplot_dplyr.sif"

prefix=$1
kraken_report=$2
alignment_report=$3
rscript=$4
taxa_names=$5

input_tsv="$prefix-rscript-input.csv"

echo -e "Species.name,num_reads,percent_classified_reads,algorithm" > $input_tsv

#echo -e "read_name,tax_name,algorithm" > $input_tsv

# trim all consecutive spaces | find all species and unclassified classifications | parses columns in correct order | sorts columns from highest percent reads to lowest
awk '{gsub(/  +/, " "); print}' $kraken_report | grep -w "[SU]" | sed -E 's|\t |\t|g' | awk -F'\t' '{print $6 "," $2 "," $1 "," "kraken"}' | sort -k3,3nr >> $input_tsv

# while IFS= read -r line; do

#     letter_code=$(echo "$line" | awk '{print $1}')
#     read_name=$(echo "$line" | awk '{print $2}')
#     tax_id=$(echo "$line" | awk '{print $3}')

#     if [[ "$letter_code" == "U" ]]; then

#         echo -e "$read_name\tunclassified\tkraken" >> $input_tsv

#     else

#         taxa_name=$(grep "^$tax_id.[|].*scientific name" $taxa_names | cut -d "|" -f 2 | tr -d '\t') # this line is grabbing the scientific name from the names.dmp taxonomy file

#         echo "$read_name,$taxa_name,kraken" >> $input_tsv

#     fi

# done < $kraken_report

#awk '{print $2 "," $3 "," "kraken"}' $kraken_report

#kraken_classified_reads=$(awk -F'\t' '$6 == "root"' $kraken_report | awk -F'\t' '{print $2}')

awk -F'\t' '{print $1 "," $2 "," $3 "," "alignment"}' $alignment_report | sort -k3,3nr >> $input_tsv

# cat $alignment_report >> $input_tsv

singularity exec --bind /scicomp $R_CONTAINER Rscript $rscript $input_tsv "$(pwd)"

mv plot.png $prefix-kraken-alignment-comparison-plot.png