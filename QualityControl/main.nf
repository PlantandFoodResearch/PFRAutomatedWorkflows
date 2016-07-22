#!/usr/bin/env nextflow


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromPath( params.global_reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.global_reads}" }  
    .set { reads }

Channel
    .fromPath( params.global_reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.global_reads}" }
    .set { reads2 }

fastqc_outdir = file(params.fastqc_outdir)


process fastqc{
    tag 'fastq'
    //executor 'lsf'
    module = params.fastqc_module
    //module 'FastQC/0.11.2'
    //cpus 8 
    //disk '10 GB'
    //memory '1000 GB'
    //time '10h'

    input:
    set file(r) from reads

    script:
        if( params.fastqc == true )
            """
            mkdir -p ${fastqc_outdir}
            fastqc -o ${fastqc_outdir}  ${r} 
            """
        else
            println "Fastqc omitted: params.fastqc!=true"
}

process fastq_screen{
    tag 'fastq_screen'
    //executor lsf
    module 'fastq_screen/v0.5.2'
    
    input:
    set file(r) from reads2
    //file cont_db

    script:
        if( params.fastq_screen == true)
            """
                fastq_screen --help
            """
        else
            error "Fastqc_screen omitted: params.fastq_screen!=true"
}

