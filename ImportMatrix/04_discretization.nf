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

    nextflow run 02_importation.nf --groupconfig 'file.conf'

    Mandatory arguments:
      --matrix                      Path to input data (must be surrounded with quotes).


    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}



// Pipeline version
version = '0.1b'
softPath = params.softPath



// Validate inputs
matrix = file(params.matrix)
if( !matrix.exists() ) exit 1, "Missing group config: '$matrix'. Specify path with --matrix"
Channel
    .fromPath( matrix )
    .ifEmpty { exit 1, "Cannot find group file\n please specify the path by using --matrix on the command line." }
    .set { matrixfile }


/*
 * STEP 1 - Divide the matrix file into multiple sub matrix
 */
process discretDivide {

    input:
    file matrixcurrent from matrixfile

    output:
    file '*.matrix' into minimatrix


    script:
    """
    $softPath/cutMatrixByLine.pl $matrixcurrent 10000000
    """
}

//to treat each sub matrix separatly
minimatrix
    .flatMap()
    .set {minimatrixsol}


/*
 * STEP 2 - Discretization
 */
process discretApply {

    input:
    file mms from minimatrixsol

    output:
    file '*.matrix' into discretizedmatrix


    script:
    prefix = mms.toString()
    """
    export OMP_NUM_THREADS=${task.cpus}
    $softPath/exportgroup.pl $mms > groupc.txt
    $softPath/indexingKmers discretizationAmeva -i $mms -o ${prefix}_discret.matrix -groups groupc.txt -buffer 1000000
    """
 }


discretizedmatrix
    .collect()
    .set {discretizedmatrixsol}

/*
 * STEP 3 - Join
 */
 process discretJoin {
  publishDir "${params.outdir}/discretization", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".matrix") > 0) "matrix/$filename"
        }

    input:
    file dms from discretizedmatrixsol

    output:
    file '*.matrix' into finalmatrix


    script:
    """
    $softPath/joinMatrixByLine.pl $dms > DISCRETmatrix.matrix
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
