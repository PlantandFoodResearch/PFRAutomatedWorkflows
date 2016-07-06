#!/usr/bin/env nextflow


/*
 * Files
 */

params.reads = "$baseDir/../KiwiTestData/*R{1,2}.fq.gz"
params.genome = "/output/genomic/plant/Actinidia/chinensis/CK51F3_01/Genome/Assembly/PseudoSanger/PS1/Pseudochromosomes/Version1/PS1.1.68/PS1.1.68.5/AllChromosomes/PS1_1.68.5.fasta"

/*
 * The reference genome file
 */

genome_file = file(params.genome)

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }


process buildGenomeIndex{
    tag 'build'
    executor 'lsf'
    module 'bowtie/1.0.0'
    /*publishDir '/output/something'*/
    cpus 1
    disk '10 GB'
    memory '1000 GB'
    time '10h'

    input:
    file genome_file

    output:
    file "genome.index*" into genome_index

    """
    bowtie-build ${genome_file} genome.index
    """
}

