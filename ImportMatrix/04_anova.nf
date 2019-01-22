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
 - 2:   discretization
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info"""
    ==========================================================================
     GECKO importation : Importation into raw matrix v${version}
    ==========================================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run 04_anova.nf --matrix 'raw.matrix' --outmatrix 'filtered.matrix' --threshold 0.05

    Mandatory arguments:
      --matrix <input file>           Path to input data (must be surrounded with quotes).
      --threshold [float]             The p-value threshold to apply ( suggested: 0.05)
    Optional:
      --outdir <output directory>       The output directory where the results will be saved

      
    """.stripIndent()
}



// Pipeline version
version = '0.1b'
softPath = params.softPath



// Validate inputs
matrix = file(params.matrix)
threshold = params.threshold


if( !matrix.exists() ) exit 1, "Missing group config: '$matrix'. Specify path with --matrix"
Channel
    .fromPath( matrix )
    .ifEmpty { exit 1, "Cannot find group file\n please specify the path by using --matrix on the command line." }
    .set { matrixfile }


/*
 * STEP 1 - Fitler the matrix file
 */
process filterMatrix {
 publishDir "${params.outdir}", mode: 'move',
        saveAs: {filename ->
            if (filename.indexOf(".matrix") > 0) "$filename"
        }
    input:
    	file matrixcurrent from matrixfile
	
	output:
		file '*.matrix' into finalmatrix

    script:
    """
    $softPath/anovaFilter $matrixcurrent filtered.matrix $threshold
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
