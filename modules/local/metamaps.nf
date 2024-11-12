process METAMAPS {
    tag "$meta.id"
    label 'process_high_memory'

    input:
    tuple val(meta), path(reads)
    path(ref)

    output:
    tuple val(meta), path("*.WIMP")                              , optional:true, emit: classification_results
    tuple val(meta), path("*species.EM.WIMP")                    , optional:true, emit: species_results
    tuple val(meta), path("*genus.EM.WIMP")                      , optional:true, emit: genus_results
    tuple val(meta), path("*reads2Taxon")                        , optional:true, emit: reads_to_taxon
    tuple val(meta), path("*krona")                              , optional:true, emit: krona_file
    tuple val(meta), path("*contigCoverage")                     , optional:true, emit: contig_coverage
    tuple val(meta), path("*lengthAndIdentitiesPerMappingUnit")  , optional:true, emit: length_and_id
    tuple val(meta), path("*final*.EM")                          , optional:true, emit: final_results

    path('*species.EM.WIMP')                                      , optional:true, emit: species_results_mqc


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/metamaps/metamaps-1.0.sif metamaps mapDirectly --all -r ${ref}/DB.fa -q ${reads} -o ./"$prefix"_classification_results --maxmemory ${params.metamaps_mem} -t ${params.metamaps_threads}
    
    singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/EMEL/singularity/metamaps/metamaps-1.0.sif metamaps classify --mappings "$prefix"_classification_results --DB ${ref} -t ${params.metamaps_threads}

    mv classification_results.EM "$prefix"_final_classification.EM

    header=\$(head -n 1 "$prefix"_classification_results.EM.WIMP)


    echo \$header > "$prefix"_classification_results_1_species.EM.WIMP

    awk '\$1 == "species"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_1_species.EM.WIMP
    
    echo \$header > "$prefix"_classification_results_2_genus.EM.WIMP

    awk '\$1 == "genus"' "$prefix"_classification_results.EM.WIMP | sort -t\$'\t' -k 4,4nr | sed 's| |_|g' | awk -v col=5 '{ \$col = \$col * 100 }1' OFS="\t" | awk -v col=6 '{ \$col = \$col * 100 }1' OFS="\t" | sed 's|_| |g' >> "$prefix"_classification_results_2_genus.EM.WIMP

    """
    
}
