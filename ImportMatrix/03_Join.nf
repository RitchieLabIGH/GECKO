#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
          GECKO MATRIX JOINTURE
========================================================================================
----------------------------------------------------------------------------------------
 
----------------------------------------------------------------------------------------
*/




def helpMessage() {
    log.info"""
    ==========================================================================
     GECKO Join Kmers : Join Kmers  v${version}
    ==========================================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run 03_Join.nf --listParams 'listParams.txt'

    Mandatory arguments:
      --listParams                  Path to input data (must be surrounded with quotes).


    Other options:
      --outdir                      The output directory where the results will be saved
    """.stripIndent()
}


listParams = file(params.listParams)
param_current = Channel.from(listParams.readLines())


// Pipeline version
version = '0.1b'

softPath = params.softPath



/*******************************************************************************
 * STEP 1 - JOIN
*******************************************************************************/
process join {
    tag "$id"
        publishDir "${params.outdir}/rawimport", mode: 'copy',
         saveAs: {filename ->
             if (filename.indexOf(".matrix") > 0) "matrix/$filename"
         }

    
        input:
        val param from param_current

        // Définit la sortie du Process
        // Tout fichier commençant par seq_ sera dans la Channel record
        output:
        file '*.matrix' into minimatrix


        script:
        """
        $softPath/JoinKMers $param
        """
}
