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
include { KRAKEN2_KRAKEN2 as KRAKEN2_STANDARD             } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_HUM                  } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_BACTERIA             } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_MAIN                 } from '../modules/nf-core/kraken2/kraken2/main'
include { FLYE                                           } from '../modules/nf-core/flye/main'
include { MEGAHIT                                        } from '../modules/nf-core/megahit/main'
include { MINIMAP2_ALIGN as ALIGN_READS                  } from '../modules/nf-core/minimap2/main'
include { MINIMAP2_ALIGN as ALIGN_READS_AGAIN            } from '../modules/nf-core/minimap2/main'
include { MINIMAP2_ALIGN as ALIGN_ASSEMBLY               } from '../modules/nf-core/minimap2/main'
include { SAMTOOLS_INDEX                                 } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT                                  } from '../modules/nf-core/samtools/sort/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS           } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { METABAT2_METABAT2                              } from '../modules/nf-core/metabat2/metabat2/main'
include { MEDAKA                                         } from '../modules/nf-core/medaka/main'
include { MAXBIN2                                        } from '../modules/nf-core/maxbin2/main'
include { BUSCO_BUSCO                                    } from '../modules/nf-core/busco/busco/main'
include { CHECKM_LINEAGEWF                               } from '../modules/nf-core/checkm/lineagewf/main'
include { BLAST_MAKEBLASTDB                              } from '../modules/nf-core/blast/makeblastdb/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'


//
// MODULE: custom, local modules
//
include { NANOPLOT as NANOPLOT_RAW                     } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_TRIMMED                 } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_TRIMMED_FILTERED        } from '../modules/local/nanoplot'
include { METAMAPS                                     } from '../modules/local/metamaps'
include { RENAME_KRAKEN_READS                          } from '../modules/local/rename_kraken_reads'
include { FASTQSCREEN                                  } from '../modules/local/fastq_screen'
include { CLASSIFIER_COMPARISON_GRAPH                  } from '../modules/local/classifier_comparison_graph'
include { SAMTOOLS_STATS                               } from '../modules/local/samtools_stats'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED    } from '../modules/local/samtools_fastq'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED      } from '../modules/local/samtools_fastq'
include { BBMAP_REFORMAT                               } from '../modules/local/bbmap_reformat'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_AGAIN       } from '../modules/local/bbmap_reformat'
include { PARSE_READS_BY_TAXON                         } from '../modules/local/parse_reads_by_taxon'
include { SPADES                                       } from '../modules/local/spades'
include { UNZIP                                        } from '../modules/local/unzip'
include { UNZIP as UNZIP_POLISHED                      } from '../modules/local/unzip'
include { ZIP                                          } from '../modules/local/zip'
include { QUAST                                        } from '../modules/local/quast'
include { QUAST as QUAST_MEDAKA                        } from '../modules/local/quast'
include { BLAST_BLASTN                                 } from '../modules/local/blastn'
include { BLAST_MEGABLAST as MEGABLAST                 } from '../modules/local/megablast'


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
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    NANOPLOT_RAW (

        INPUT_CHECK.out.reads

    )

    PORECHOP_PORECHOP (

        INPUT_CHECK.out.reads    

    )

    CHOPPER (

        PORECHOP_PORECHOP.out.reads

    )

    NANOPLOT_TRIMMED (

        CHOPPER.out.fastq

    )

    trimmed_reads = CHOPPER.out.fastq

    // host/human read removal with minimap2
    if (params.alignment_based_filtering) {

        ALIGN_READS (

            trimmed_reads,
            [[params.minimap2_meta],[params.minimap2_ref]],
            true,
            false,
            false        

        )

        SAMTOOLS_STATS (

            ALIGN_READS.out.bam

        )

        // capturing all unaligned reads and converting back into a fastq
        SAMTOOLS_FASTQ_UNMAPPED (

            ALIGN_READS.out.bam,
            false

        )

        // using bbmap suite to remove empty reads and deduplicate reads
        BBMAP_REFORMAT (

            SAMTOOLS_FASTQ_UNMAPPED.out.fastq

        )

        // capturing aligned reads (human)
        SAMTOOLS_FASTQ_MAPPED (

            ALIGN_READS.out.bam,
            false

        )

        BBMAP_REFORMAT_AGAIN (

            SAMTOOLS_FASTQ_MAPPED.out.fastq

        )

        // setting filtered_reads channel equal to unaligned reads
        filtered_reads = BBMAP_REFORMAT.out.fastq

        // running nanoplot again to compare read stat pre and post filter
        NANOPLOT_TRIMMED_FILTERED (

            filtered_reads

        )

    } else {

        filtered_reads = CHOPPER.out.fastq

    }

    // taxonmic profiling with metamaps
    if (params.use_metamaps) {

        METAMAPS (
            filtered_reads,
            params.metamaps_db
        )

    }

    // taxonmic profiling with kraken2
    if (params.use_kraken2) {
        
        KRAKEN2_MAIN (

            filtered_reads,
            params.kraken_db_main,
            true,
            false

        )

    }

    // contaminant screenin and limited taxonmic profiling with fastqscreen
    FASTQSCREEN (
        filtered_reads,
        params.fastq_screen_conf
    )

    // horizontal bar chart comparison of metamaps and kraken2 results
    if (params.use_kraken2 && params.use_metamaps) {
        
        graph_input = METAMAPS.out.species_results.join(KRAKEN2_MAIN.out.report)

        CLASSIFIER_COMPARISON_GRAPH (
            graph_input
        )

    }

    // parse out reads classified under a designated tax ID for further analyses
    if (params.parse_reads_by_taxon && params.use_metamaps) {

        prbt_input = METAMAPS.out.reads_to_taxon.join(filtered_reads)

        PARSE_READS_BY_TAXON (

            prbt_input

        )

    }
    
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


    // if (params.use_minimap2) {


    //     ALIGN_ASSEMBLY (

    //         filtered_reads,
    //         unzip_channel,
    //         true,
    //         false,
    //         false

    //     )

    //     // samtools steps below are only needed if using metabat2 for binning
    //     // SAMTOOLS_SORT (

    //     //     ALIGN_ASSEMBLY.out.bam,
    //     //     unzip_channel

    //     // )

    //     // SAMTOOLS_INDEX (

    //     //     SAMTOOLS_SORT.out.bam

    //     // )


    // }


    //assembly qc with quast
    QUAST (

        unzip_channel, // consensus (one or more assemblies)

    )


    if (params.use_medaka) {

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

        // metabat2_sum_input_ch = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)

        // METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS (

        //     metabat2_sum_input_ch

        // )

        // metabat2_metabat2_input_ch = unzip_channel.join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth)

        // METABAT2_METABAT2 (

        //     metabat2_metabat2_input_ch

        // )


        // CHECKM_LINEAGEWF (

        //     MAXBIN2.out.binned_fastas,
        //     "fasta.gz",
        //     []

        // )

        BUSCO_BUSCO (

            MAXBIN2.out.binned_fastas,
            "genome",
            "eukaryota_odb10",
            [],
            []

        )


    }

    if (!params.skip_blast) {

        BLAST_BLASTN (
            unzip_channel,
            params.blast_db
        )

        // MEGABLAST (

        //     unzip_channel,
        //     params.blast_standard_db

        // )

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

    //ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_RAW.out.txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_TRIMMED.out.txt.collect{it[1]}.ifEmpty([]))
    if (params.use_kraken2) {
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_MAIN.out.report_mqc.collect().ifEmpty([]))
    }
    
    if (params.use_medaka) {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST_MEDAKA.out.report.collect().ifEmpty([]))
    } else {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))
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
