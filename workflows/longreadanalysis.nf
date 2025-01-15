/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

//skip this step to prevent "--fasta not specified" errors
//WorkflowLongreadanalysis.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                         } from '../modules/nf-core/fastqc/main'
include { PORECHOP_PORECHOP                              } from '../modules/nf-core/porechop/porechop/main'
include { CHOPPER                                        } from '../modules/nf-core/chopper/main'
include { NONPAREIL_NONPAREIL                            } from '../modules/nf-core/nonpareil/nonpareil/main'
include { NONPAREIL_CURVE                                } from '../modules/nf-core/nonpareil/curve/main'
include { NONPAREIL_NONPAREILCURVESR                     } from '../modules/nf-core/nonpareil/nonpareilcurvesr/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_MAIN                } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_PROTEIN             } from '../modules/nf-core/kraken2/kraken2/main'
include { KRONA_KRONADB                                  } from '../modules/nf-core/krona/krona_db/main'
include { KRAKENTOOLS_KREPORT2KRONA                      } from '../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS                 } from '../modules/nf-core/krakentools/extractkrakenreads/main'
include { KRONA_KTIMPORTTEXT as KRONA_KRAKEN             } from '../modules/nf-core/krona/ktimporttext/main'
include { KRONA_KTIMPORTTEXT as KRONA_METAMAPS           } from '../modules/nf-core/krona/ktimporttext/main'
include { FLYE                                           } from '../modules/nf-core/flye/main'
include { MEGAHIT                                        } from '../modules/nf-core/megahit/main'
include { MINIMAP2_ALIGN as ALIGN_READS                  } from '../modules/nf-core/minimap2/main'
include { SAMTOOLS_DEPTH                                 } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_INDEX                                 } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COVERAGE                              } from '../modules/nf-core/samtools/coverage/main'
include { MEDAKA                                         } from '../modules/nf-core/medaka/main'
include { MAXBIN2                                        } from '../modules/nf-core/maxbin2/main'
include { CHECKM2_PREDICT                                } from '../modules/nf-core/checkm2/predict/main'
include { BUSCO_BUSCO                                    } from '../modules/nf-core/busco/busco/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'


//
// MODULE: custom, local modules
//
include { NANOPLOT as NANOPLOT_RAW                        } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_TRIMMED                    } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_KRAKEN_TAXON_FILTERED      } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_ALIGNMENT_TAXON_FILTERED   } from '../modules/local/nanoplot'
include { METAMAPS                                        } from '../modules/local/metamaps'
include { FASTQSCREEN                                     } from '../modules/local/fastq_screen'
include { CLASSIFIER_COMPARISON_GRAPH                     } from '../modules/local/classifier_comparison_graph'
include { SAMTOOLS_STATS                                  } from '../modules/local/samtools_stats'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED       } from '../modules/local/samtools_fastq'
include { ALIGNMENT_CLASSIFY                              } from '../modules/local/alignment_classify'
include { SAMTOOLS_VIEW as SAMTOOLS_QUALITY_FILTER        } from '../modules/local/samtools_view'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         } from '../modules/local/samtools_fastq'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_SUBSAMPLE      } from '../modules/local/bbmap_reformat_subsample'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_CLEAN_UNMAPPED } from '../modules/local/bbmap_reformat'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_CLEAN_MAPPED   } from '../modules/local/bbmap_reformat'
include { PARSE_READS_BY_TAXON                            } from '../modules/local/parse_reads_by_taxon'
include { SPADES                                          } from '../modules/local/spades'
include { UNZIP                                           } from '../modules/local/unzip'
include { UNZIP as UNZIP_POLISHED                         } from '../modules/local/unzip'
include { ZIP                                             } from '../modules/local/zip'
include { QUAST                                           } from '../modules/local/quast'
include { QUAST as QUAST_MEDAKA                           } from '../modules/local/quast'
include { BLAST_BLASTN                                    } from '../modules/local/blastn'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

ch_versions = Channel.empty()
ch_multiqc_files = Channel.empty()

workflow LONGREADANALYSIS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    if (!params.skip_subsample) {

        BBMAP_REFORMAT_SUBSAMPLE (
            INPUT_CHECK.out.reads,
            params.num_subsamples
        )

        raw_reads = BBMAP_REFORMAT_SUBSAMPLE.out.fastq

    } else {

        raw_reads = INPUT_CHECK.out.reads

    }

    NANOPLOT_RAW (

        raw_reads

    )

    if (params.skip_porechop) {

        trimmed_reads = raw_reads
    
    } else {

        PORECHOP_PORECHOP (

            raw_reads    

        )

        trimmed_reads = PORECHOP_PORECHOP.out.reads

    }

    CHOPPER (

        trimmed_reads

    )

    NANOPLOT_TRIMMED (

        CHOPPER.out.fastq

    )

    trimmed_reads = CHOPPER.out.fastq


    if (!params.skip_nonpareil) {
        
        NONPAREIL_NONPAREIL (

            trimmed_reads,
            "fastq",
            "alignment"

        )

        NONPAREIL_CURVE (

            NONPAREIL_NONPAREIL.out.npo

        )

        // npo_reports = Channel.empty()
        // npo_reports = npo_reports.mix(NONPAREIL_NONPAREIL.out.npo.collect().ifEmpty([]))
        // npo_reports_tuple = npo_reports
        //                     .map {meta, npo -> [[id: 'all'], npo] }
        //                     .groupTuple()

        NONPAREIL_NONPAREILCURVESR (

            NONPAREIL_NONPAREIL.out.npo,
            NONPAREIL_CURVE.out.png

        )

    }
    

    // host/human read removal with minimap2
    if (!params.skip_alignment_based_filtering) {

        ALIGN_READS (

            trimmed_reads,
            [[params.minimap2_meta],[params.minimap2_index]],
            true,
            false,
            false        

        )

        // capturing all unaligned reads and converting back into a fastq
        SAMTOOLS_FASTQ_UNMAPPED (

            ALIGN_READS.out.bam,
            false

        )

        // using bbmap suite to remove empty reads and deduplicate reads
        BBMAP_REFORMAT_CLEAN_UNMAPPED (

            SAMTOOLS_FASTQ_UNMAPPED.out.fastq

        )

        // custom module that parses a bam file and pulls out each alignment that corresponds to a reference sequence under a specified tax id
        ALIGNMENT_CLASSIFY (

            ALIGN_READS.out.bam,
            params.seq2tax_map,
            params.my_tax_ids,
            params.include_children

        )

        SAMTOOLS_STATS (

            ALIGNMENT_CLASSIFY.out.bam

        )

        SAMTOOLS_DEPTH (

            ALIGNMENT_CLASSIFY.out.bam

        )

        SAMTOOLS_INDEX (

            ALIGNMENT_CLASSIFY.out.bam

        )

        SAMTOOLS_COVERAGE (

            ALIGNMENT_CLASSIFY.out.bam.join(SAMTOOLS_INDEX.out.bai)

        )

        // filtering for reads with a mapping quality at or above params.mapping_quality
        SAMTOOLS_QUALITY_FILTER (

            ALIGNMENT_CLASSIFY.out.bam,
            false

        )

        // capturing aligned reads
        SAMTOOLS_FASTQ_MAPPED (

            SAMTOOLS_QUALITY_FILTER.out.bam,
            false

        )

        // using bbmap suite to remove empty reads and deduplicate aligned reads
        BBMAP_REFORMAT_CLEAN_MAPPED (

            SAMTOOLS_FASTQ_MAPPED.out.fastq

        )

        // setting filtered_reads channel equal to aligned reads
        filtered_reads = BBMAP_REFORMAT_CLEAN_MAPPED.out.fastq

        // running nanoplot again to compare read stat pre and post filter
        NANOPLOT_ALIGNMENT_TAXON_FILTERED (

            filtered_reads

        )

    } else {

        filtered_reads = CHOPPER.out.fastq

    }

    // taxonmic profiling with metamaps
    if (!params.skip_metamaps) {

        METAMAPS (
            CHOPPER.out.fastq,
            params.metamaps_db
        )

    }

    // taxonmic profiling with kraken2
    if (!params.skip_kraken2) {
        
        KRAKEN2_MAIN (

            CHOPPER.out.fastq,
            params.kraken_db_main,
            true,
            true

        )

        KRAKENTOOLS_KREPORT2KRONA (
            KRAKEN2_MAIN.out.report
        )

        KRONA_KRAKEN (
            KRAKENTOOLS_KREPORT2KRONA.out.txt
        )

    }

    if (!params.skip_fastq_screen) {

        // contaminant screenin and limited taxonmic profiling with fastqscreen
        FASTQSCREEN (
            CHOPPER.out.fastq,
            params.fastq_screen_conf
        )

    }

    // horizontal bar chart comparison of metamaps and kraken2 results
    if (!params.skip_kraken2 && !params.skip_metamaps) {
        
        graph_input = METAMAPS.out.species_results.join(KRAKEN2_MAIN.out.report)

        CLASSIFIER_COMPARISON_GRAPH (
            graph_input
        )

    }

    // parse out reads classified under a designated tax ID
    if (!params.skip_kraken2_parse_reads_by_taxon && !params.skip_kraken2) {

        KRAKENTOOLS_EXTRACTKRAKENREADS (

            params.kraken_tax_ids,
            KRAKEN2_MAIN.out.classified_reads_assignment,
            KRAKEN2_MAIN.out.classified_reads_fastq,
            KRAKEN2_MAIN.out.report

        )

        NANOPLOT_KRAKEN_TAXON_FILTERED (

            KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads

        )

        filtered_reads = KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads

    }

    if (!params.skip_kraken2_protein) {
    
        KRAKEN2_PROTEIN (

            filtered_reads,
            params.kraken_db_protein,
            true,
            false

        )

    }

    if (!params.skip_filter_by_kraken2_protein && !params.skip_kraken2_protein) {

        filtered_reads = KRAKEN2_PROTEIN.out.classified_reads_fastq

    }

    if (!params.skip_assembly) {
    
        if (params.assembler == "spades") {
            
            spades_map = filtered_reads.map { meta, fastq -> [ meta, [], [], fastq ] }

            contigs_produced = true

            SPADES (

                spades_map,
                [],
                [],         

            )

            assembly_ch = SPADES.out.contigs

        } else if (params.assembler == "flye") {

            FLYE (

                filtered_reads,
                '--nano-hq'

            )

            assembly_ch = FLYE.out.fasta

        } else if (params.assembler == "megahit") {

            MEGAHIT (

                filtered_reads

            )

            assembly_ch = MEGAHIT.out.contigs
        
        }

        UNZIP (

            assembly_ch

        )

        unzip_channel = UNZIP.out.unzip_contigs


        //assembly qc with quast
        QUAST (

            unzip_channel, // consensus (one or more assemblies)

        )


        if (!params.skip_medaka) {

            medaka_input = filtered_reads.join(unzip_channel)

            //polishing contigs with medaka
            MEDAKA (

                medaka_input

            )

            assembly_ch = MEDAKA.out.assembly

            QUAST_MEDAKA (

                assembly_ch

            )
        
            //generated unzipped version of the spades contigs
            UNZIP_POLISHED (

                assembly_ch

            )

            unzip_channel = UNZIP_POLISHED.out.unzip_contigs

        }

        if (!params.skip_binning) {

            maxbin_map = unzip_channel.join(filtered_reads)

            MAXBIN2 (

                maxbin_map

            )

            if (params.bin_QC == "checkm") {

                CHECKM2_PREDICT (

                    MAXBIN2.out.binned_fastas,
                    [["checkm2-DB"],["${params.checkm2_db}"]]

                )

            } else {

                BUSCO_BUSCO (

                    MAXBIN2.out.binned_fastas,
                    "genome",
                    "bacteria_odb10",
                    [],
                    []

                )
            
            }


        }

        if (!params.skip_blast) {

            BLAST_BLASTN (
                unzip_channel,
                params.blast_db
            )

        }
    }
    
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )




    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowLongreadanalysis.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowLongreadanalysis.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_RAW.out.txt.collect{it[1]}.ifEmpty([]))
    if (!params.skip_porechop) {
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log.collect().ifEmpty([]))
    }
    
    if (!params.skip_nonpareil) {
        ch_multiqc_files = ch_multiqc_files.mix(NONPAREIL_NONPAREILCURVESR.out.mqc_json.collect().ifEmpty([]))
    }

    if (!params.skip_alignment_based_filtering) {
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.txt.map{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_DEPTH.out.tsv.map{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_COVERAGE.out.coverage.map{it[1]})
    }

    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_TRIMMED.out.txt.collect{it[1]}.ifEmpty([]))
    if (!params.skip_kraken2) {
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_MAIN.out.report_mqc.collect().ifEmpty([]))
    }
    
    if (!params.skip_medaka) {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST_MEDAKA.out.report.collect().ifEmpty([]))
    } else if (!params.skip_assembly) {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))
    }

    if (!params.skip_binning && params.bin_QC != "checkm" ) {
        ch_multiqc_files = ch_multiqc_files.mix(BUSCO_BUSCO.out.short_summaries_multiqc.collect().ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
