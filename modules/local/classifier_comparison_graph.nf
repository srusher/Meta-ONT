process CLASSIFIER_COMPARISON_GRAPH {
    tag "$meta.id"

    input:
    tuple val(meta), path(mmap_report), path(k2_report)

    output:
    path('*.png')         , optional:true, emit: comparison_plot

    script:
    def prefix = "${meta.id}"
    """
    bash ${projectDir}/bin/process_metamaps_results.sh $mmap_report > ${prefix}_mmap_data
    bash ${projectDir}/bin/process_kraken2_results.sh $k2_report > ${prefix}_k2_data

    sort -t\$'\t' -k 2,2nr ${prefix}_mmap_data | awk 'NR <= 15' > ${prefix}_merged_sorted
    sort -t\$'\t' -k 2,2nr ${prefix}_k2_data | awk 'NR <= 15' >> ${prefix}_merged_sorted

    sort -t\$'\t' -k 2,2nr ${prefix}_merged_sorted > ${prefix}_temp
    cat ${prefix}_temp > ${prefix}_merged_sorted

    echo -e "Species\tAbundance(%)\tClassifier" > ${prefix}_merged_sorted_wheader
    cat ${prefix}_merged_sorted >> ${prefix}_merged_sorted_wheader

    singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/seaborn/seaborn-latest.sif python3 ${projectDir}/bin/meta_graph.py "${prefix}_merged_sorted_wheader" "${prefix}"
    """
}
