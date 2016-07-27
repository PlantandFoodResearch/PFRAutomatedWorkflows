#!/usr/bin/env nextflow


/*
 * Files
 */

/*
 * The reference genome file
 */

genome = Channel.fromPath(params.genome)


/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }


process testData{
    tag 'test'
    executor 'lsf'
    
    input:
    set pair_id, file(reads) from read_pairs

    script:
    template 'test_data.sh'
}


process buildGenomeIndex{
    tag 'build'
    executor 'lsf'
    module 'bowtie/1.0.0'
    publishDir 'output'
    cpus 1
    disk '10 GB'
    memory '1000 GB'
    time '10h'

    input:
    file species from genome

    output:
    file "${dbName}*" into genome_index

    script:
    dbName = species.baseName

    """
    bowtie-build ${species} genome_index
    """
}




