#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          GECKO MATRIX
========================================================================================
----------------------------------------------------------------------------------------
 
 
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
      --matrixDiscrete                      Path to input data (must be surrounded with quotes).
      --matrixReal                      Path to input data (must be surrounded with quotes).


    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}


// Pipeline version
version = '0.1b'
softPath = params.softPath



// Validate inputs
matrixDiscrete = file(params.matrixDiscrete)
matrixReal     = file(params.matrixReal)


if( !matrixDiscrete.exists() ) exit 1, "Missing group config: '$matrixDiscrete'. Specify path with --matrixDiscrete"
if( !matrixReal.exists() ) exit 1, "Missing group config: '$matrixReal'. Specify path with --matrixReal"
Channel
    .fromPath( matrixDiscrete )
    .ifEmpty { exit 1, "Cannot find group file\n please specify the path by using --matrixDiscrete on the command line." }
    .set { matrixDiscretefile }
Channel
    .fromPath( matrixReal )
    .ifEmpty { exit 1, "Cannot find group file\n please specify the path by using --matrixReal on the command line." }
    .set { matrixRealfile }





process toreal {
    publishDir "${params.outdir}/filtering", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".matrix") > 0) "final/$filename"
        }
    
    input:
    file mdis from matrixDiscretefile
    file mreal from matrixRealfile

    output:
    file '*.matrix' into outpumatrix


    script:
    """
    $softPath/reduct2realFile.pl $mdis $mreal > FILTEREDmatrix_RealCounts.matrix
    """
}


process toMLformat {
    publishDir "${params.outdir}/filtering", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".matrix") > 0) "final/$filename"
        }
    
    input:
    file mreal from outpumatrix

    output:
    file '*.matrix' into outpumatrix


    script:
    """
    $softPath/transformIntoML $mreal > FILTEREDmatrix_RealCounts_MLformat.matrix
    """
}
