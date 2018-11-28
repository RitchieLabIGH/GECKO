#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          GECKO MATRIX
========================================================================================
----------------------------------------------------------------------------------------
 Pipeline overview:
 - 1:   cut configuration file into pieces
 - 2:   import pieces
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info"""
    ==========================================================================
     GECKO importation : Importation into raw matrix v${version}
    ==========================================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run 02_importation.nf --groupconfig 'file.conf'

    Mandatory arguments:
      --groupconfig                       Path to input data (must be surrounded with quotes).


    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1b'

softPath = params.softPath

// Validate inputs
groupconfig = file(params.groupconfig)
if( !groupconfig.exists() ) exit 1, "Missing group config: '$groupconfig'. Specify path with --groupconfig"
Channel
    .fromPath( groupconfig )
    .ifEmpty { exit 1, "Cannot find group file\n please specify the path by using --groupconfig on the command line." }
    .first()
    .set { groupconfigfile }


/*
 * STEP 1 - Divide the configuration file into multiple files
 */
process rawsubextract {
    publishDir "${params.outdir}/rawimport", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".conf") > 0) "conf/$filename"
        }

    input:
    file groupfile from groupconfigfile

    output:
    file '*.conf' into miniconf


    script:
    """
    split -l 5 --additional-suffix=.conf $groupfile subgroup
    """
}

miniconf
    .flatMap()
    .set {miniconfsol}


/*
 * STEP 2 - Import sub matrix
 */
 process rawsubimport {
     publishDir "${params.outdir}/rawimport", mode: 'copy',
         saveAs: {filename ->
             if (filename.indexOf(".conf") > 0) "matrix/$filename"
         }

     input:
     file groupfile from miniconfsol

     output:
     file '*.matrix' into minimatrix


     script:
     prefix = groupfile.toString()
     """
     $softPath/indexingKmers importation -i $groupfile -o ${prefix}.matrix -threshold 5
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
