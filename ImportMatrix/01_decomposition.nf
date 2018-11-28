#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          GECKO MATRIX
========================================================================================
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   FastQC for raw sequencing reads quality control
 - 2:   Trim Galore! for adapter trimming
 - 3.1: k-mer
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info"""
    ==========================================================================
     GECKO decomposition : Decomposition into k-mers v${version}
    ==========================================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run 01_decomposition.nf --reads '*_R{1,2}.fastq.gz'

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).

    Options:
      --singleEnd                   Specifies that the input is single end reads

    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1b'

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0


/*
 * Create a channel for input read files
 */
params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { raw_reads_fastqc; raw_reads_trimgalore }


softPath = params.softPath



/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    """
}



/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
    trimgalore_fastqc_reports = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from raw_reads_trimgalore

        output:
        file '*.fq' into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            $softPath/trim_galore --dont_gzip --fastqc $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            $softPath/trim_galore --dont_gzip --paired --fastqc $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}




/*
 * STEP 3 - Jellyfish count
 */
process jellyfish {
    tag "$prefix"
    publishDir "${params.outdir}/jellyfish", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("ojf") > 0) "binary/$filename"
        }

    input:
    file reads from trimmed_reads

    output:
    file '*.ojf' into jellyfish_ojf

    script:
    prefix = reads[0].toString() - ~/(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    jellyfish count -C -m 30 -p 256 --disk -t ${task.cpus} -c 64 -s 500000000 -o ${prefix}.ojf $reads
    """
}



/*
 * STEP 4 - Jellyfish extract
 */
process extract {
    tag "$prefix"
    publishDir "${params.outdir}/jellyfish", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("ojf.tab") > 0) "text/$filename"
        }

    input:
    file jellyin from jellyfish_ojf

    output:
    file '*.ojf.tab' into jellyfish_tab

    script:
    prefix = jellyin.toString()
    """
    jellyfish dump -c -t -o ${prefix}.tab $jellyin
    """
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution project held at ${workflow.projectDir} "
    println "Execution launched at ${workflow.launchDir} "
    """
    echo "done done done"
    """
}
