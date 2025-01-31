#!/bin/bash

# POC: Sam Rusher (rtq0@cdc.gov)

# DESCRIPTION: This script is designed perform basecalling on raw ONT reads using Dorado. It can take in fast5 or pod5 files as inputs. If the file is in fast5 format the script will convert it into pod5 before running through dorado

#$ -m abe 
#$ -o qsub_logs/ 
#$ -e qsub_logs/
#$ -M rtq0@cdc.gov
#$ -N Dorado
#$ -l "gpu=1 h_vmem=32G" 
#$ -q gpu.q
#$ -cwd

for file in "$@"; do

    if [[ $file == *".gz" ]]; then

        echo -e "\nInput file must be unzipped/decompressed\n"
        break

    fi

    if [[ -f $file ]]; then

        if [[ "$file" == *".fast5" ]]; then

            output_file="$(basename "$file" .fast5)"_converted.pod5

            echo -e "\nConverting fast5 into pod5...\n"

            singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/pod5/pod5.sif pod5 convert fast5 $file $output_file
        
        else 

            output_file=$file

        fi

        bam_file="$(basename "$output_file" .pod5)".bam

        barcode_output="$(echo $output_file | cut  -d '.' -f1)"

        barcode_output="$(basename $barcode_output)"

        echo -e "\nBasecalling with dorado...\n"

        singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/dorado/dorado-0.7.3_build.sif dorado basecaller /scicomp/home-pure/rtq0/EMEL-GWA/singularity/dorado/basecalling-models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 $output_file --kit-name SQK-RPB114-24 --no-trim --device cuda:all > $bam_file

        # demultiplexing basecalled reads so that they are sorted by their original barcode label
        singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/dorado/dorado-0.7.3_build.sif dorado demux --no-classify --output-dir ./by-barcode/$barcode_output $bam_file

        mkdir -p ./by-barcode/$barcode_output/trimmed

        echo -e "\nBasecalling complete\n"

        barcodes="$(find ./by-barcode/$barcode_output -name "*.bam" | rev | cut -d "_" -f1 | rev | cut -d 'e' -f2 | cut -d '.' -f1 | awk '!/d/' | sort -n | uniq)"
        barcodes="$barcodes unclassified"

        for i in $barcodes; do

            if [[ $i == "unclassified" ]]; then

                bams="$(find ./by-barcode/$barcode_output -name "*$i.bam")"

            else

                bams="$(find ./by-barcode/$barcode_output -name "*barcode$i.bam")"
            
            fi

            for j in $bams; do

                trim_output="$(basename $j)"

                # trimming adapters and primers AFTER demultiplex step - otherwise trimming could interfere with barcode classification
                singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/dorado/dorado-0.7.3_build.sif dorado trim $j > ./by-barcode/$barcode_output/trimmed/$trim_output

            done

        done


        barcodes="$(find ./by-barcode/$barcode_output/trimmed -name "*.bam" | rev | cut -d "_" -f1 | rev | cut -d 'e' -f2 | cut -d '.' -f1 | awk '!/d/' | sort -n | uniq)"
        barcodes="$barcodes unclassified"

        if [[ ! -d ./fastqs ]]; then

            mkdir ./fastqs

        fi

        for i in $barcodes; do

            if [[ $i == "unclassified" ]]; then

                bams="$(find ./by-barcode/$barcode_output/trimmed -name "*$i.bam")"
                output_dir=./fastqs/$i.fastq

            else

                bams="$(find ./by-barcode/$barcode_output/trimmed -name "*barcode$i.bam")"
                output_dir=./fastqs/barcode$i.fastq
            
            fi

            for j in $bams; do

                singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/samtools/samtools:1.9--h91753b0_8 samtools fastq $j >> $output_dir

            done

        done

        # fastq_file="$(basename "$bam_file" .bam)".fastq

        # singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/samtools/samtools:1.9--h91753b0_8 samtools fastq $bam_file > $fastq_file

        #rm -f $bam_file

    else

        echo -e "\nThe file: $file, DOES NOT EXIST -- SKIPPING\n"

    fi

done