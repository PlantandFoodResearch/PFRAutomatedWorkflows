/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'RNA-Toy'.
 *
 *   RNA-Toy is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNA-Toy is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNA-Toy.  If not, see <http://www.gnu.org/licenses/>.
 */
 
/* 
 * Proof of concept Nextflow based RNA-Seq pipeline
 * 
 * Authors:
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Emilio Palumbo <emiliopalumbo@gmail.com> 
 */ 

 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */

/*
params.reads ="/workspace/cfpcxs/PFRAutomatedWorkflows/KiwiTestData/*.R{1,2}.fq"
params.annot = "/workspace/cfpcxs/KiwiTestData/kiwitest.gff3"
params.genome = "/workspace/cfpcxs/PFRAutomatedWorkflows/KiwiTestData/kiwitest.fasta"
params.outdir = 'results13'
*/

params.reads ="/workspace/cfpcxs/PFRAutomatedWorkflows/RNAseq/Petunia/*_R{1,2}.fastq" /*Defines the path for the reads file */
params.annot = "/workspace/cfpcxs/PFRAutomatedWorkflows/RNAseq/Petunia/annotation.gff" /* Defines the path for the annotation file */
params.genome = "/workspace/cfpcxs/PFRAutomatedWorkflows/RNAseq/Petunia/genome.fasta"  /* Defines the part for the reference genome file (fasta format) */
params.outdir = 'results16'   /* Defines an output directory where the you will be able to see the output files/data of this pipeline */

log.info """\
         R N A T O Y   P I P E L I N E    
         =============================
         genome: ${params.genome}
         annot : ${params.annot}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()

/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)

/*The reference genome file for contamination*/
/* It will be needed/needing for sort the reads and remove the contaminated reads */

SMRNA="/software/bioinformatics/sortmerna-2.1-linux-64"
DB="${SMRNA}/rRNA_databases"
INDEX="${SMRNA}/index"

/* Put hte pathway of each database contaminatiom genome (fasta fornat) with their index (db format): if no index, you need to create them with the SortMeRNA command */

SORTMERNADB="\
${DB}/silva-bac-16s-id90.fasta,${INDEX}/silva-bac-16s-db:\
${DB}/silva-bac-23s-id98.fasta,${INDEX}/silva-bac-23s-db:\
${DB}/silva-arc-16s-id95.fasta,${INDEX}/silva-arc-16s-db:\
${DB}/silva-arc-23s-id98.fasta,${INDEX}/silva-arc-23s-db:\
${DB}/silva-euk-18s-id95.fasta,${INDEX}/silva-euk-18s-db:\
${DB}/rfam-5s-database-id98.fasta,${INDEX}/rfam-5s-db:\
${DB}/rfam-5.8s-database-id98.fasta,${INDEX}/rfam-5.8s-db"
 
/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */

/*Def PairChannel create the channel with the tuples above, it allows to use twice or more the channel in different process */

def PairChannel(parameterR){
Channel
    .fromFilePairs( parameterR )
    .ifEmpty { error "Cannot find any reads matching: ${parameterR}" }
}

/*Test on small sample of real data: test if the pipeline works*/
 
process sample {
       module "seqtk"
       publishDir params.outdir
       executor "lsf"

        input:
        set pair_id, file(reads) from PairChannel(params.reads)

        output:
        set pair_id, file("*.sub.fq") into Sample

"""
seqtk sample -s100 ${pair_id}_R1.fastq 10000 > ${pair_id}.R1.sub.fq
seqtk sample -s100 ${pair_id}_R2.fastq 10000 > ${pair_id}.R2.sub.fq
"""
}


/* Test command Fastqc: take Fastq files which contains the reads and return Quality control report of these reads in html and zip format. There is the first check of the quality of data juste after the sequencing. */

process qualitycontrol {
        module "FastQC"
	publishDir params.outdir, mode: 'copy'   /* Copy the return data in the output directory: allow to see the data files  */
	executor "lsf" 

        input:
        set pair_id, file(reads) from Sample /*Sample contains small samples of the reads data: the pair_id is an ID in common for each forward and reverse read, for example if you have Kiwitestdata1.R1.fq, Kiwitestdata1.R2.fq, the pair_id will be "Kiwitestdata1", file(reads) just take the corresponding reads of the Sample channel */

        output:
        set "*.zip", "*.html" into qualitycontrolR /* Get back the QCreports returned by the fastqc command in the channel QualitycontrolR */
        set pair_id, file(reads) into Sample2  /* In Nextflow, you can't use the same channel twice, so as we need to use the reads containing in the channel sample in other process, we need to stock the reads data in other channel, for example: Sample2 */
"""
fastqc ${reads}  /* fastqc is a command which take the reads and retunr QCreports*/

"""
}

/* Test TRIMMING */
/* In the next steps, we will check the quality of reads and if necessary make some trimming to enhance the quality, before to go on the alignment step */

/* Test command SortMeRNA: it allows to remove the contaminated reads. We need to merge the reads before using the SortMeRNA filter, with the bash script "merge-paired-reads" we obtain merged/interleaved fastq files. With these files we make the SortMeRNA and getabck only the data which are not contaminated. We need to unmerge the data with the script "unmerge-paired-reads" to obtain fastq files unmerged in the aim to use them in next steps*/

process SortMeRNA {
	publishDir params.outdir, mode: 'copy'
	module "sortmerna/2.1"
	executor "lsf"

	input:
	set pair_id, file(reads) from Sample2

	output:
	file('*.MERGED.fastq') into MergedData
	file('*_filtered.fastq') into FINTERLEAVEDFiltered
	set pair_id, file('*_R*.fq') into FQFiltered  /* Here, we get back the final unmerged data sorted in the FQFiltered channel: reads file in fastq format without the contaminated data */


/* Probably need to gunzip with zcat <(zcat filename.gz) */

"""
bash merge-paired-reads.sh ${pair_id}.R1.sub.fq ${pair_id}.R2.sub.fq ${pair_id}.MERGED.fastq /*Merged the forward and reverse read in one file*/
sortmerna --ref ${SORTMERNADB} --reads ${pair_id}.MERGED.fastq --paired_in --fastx --aligned ${pair_id}_rRNA --other ${pair_id}_filtered /* Aligned the reads data with different contamination genome inquired in the "--ref", the variable "{SORTMERNADB}" need to be inquired by the user on the top of the pipeline by the user. The aligned data are the contaminated data so we remove them and just get back the data which were not aligned with the contamination genome so the data ouput in " --other" */ 
bash unmerge-paired-reads.sh ${pair_id}_filtered.fastq ${pair_id}_R1.fq ${pair_id}_R2.fq /*Unmerged the sorted data*/ 
"""
}



/*Test command Fastqc on the filtered fastq*/
/* We need to check if the quality of data will not decrease after the sorting, and we can continue to use it in other steps */

process QualityCFiltered {
        module "FastQC"
        publishDir params.outdir, mode: 'copy'
	executor "lsf"

        input:
        set pair_id, file(freads) from FQFiltered

        output:
        set "*.zip", "*.html" into qualitycontrolRF
        set pair_id, file(freads) into FQFiltered2
"""
fastqc ${freads}

"""
} 


/*Test command Trimmomatic*/
/* Sometimes we need to make some trimming in the data for inhance the quality of data, here we use the trimmomatic script */

process Trimmomatic {
        module "java"
        module "Trimmomatic-0.36"
        publishDir params.outdir, mode: 'copy'
	executor "lsf"

        input:
        set pair_id, file(fqreads) from FQFiltered2

        output:
        set pair_id, file('*trimmomatic_*P.fq.gz') into TrimmomaticPaired /* We get back the data which are correctly trimmed (have survived of the trimming): both of the reads (forward and reverse) have survived of the trimming, ine the TrimmomaticPaired channel */
        set pair_id, file('*trimmomatic_*U.fq.gz') into TrimmomaticUnpaired /*  We get back the data which are not correctly trimmed: at least one of the paired-reads (forward or reverse) or the two have be removed by the trimming, ine the TrimmomaticUnPaired channel *

"""
java -jar /software/bioinformatics/Trimmomatic-0.36/trimmomatic.jar PE ${pair_id}_R1.fq  ${pair_id}_R2.fq -baseout ${pair_id}.trimmomatic.fq.gz ILLUMINACLIP:"/software/bioinformatics/Trimmomatic-0.36/adapters/TruSeq2-PE.fa":2:30:10 SLIDINGWINDOW:5:20 MINLEN:50
/* We need to indicate if it's pair-ended reads with "PE", just after that we put the reads data, in the -baseout it will return the ouput file for each pair ended read so the Paired and UnPaired Trimmomatic for forward and reverse reads so 2 file for Paired (R1Paired, R2Paired) and 2 files for UnPaired (R1UnPaired, R2UnPaired) */
/* Illuminaclip: find and remove Illumina adapters which could contaminated the reads data during the sequencing step: it's need a fasta file containing illumina adapters (here we inquire the pathway of this file), the "30" is for clipping when the score is lower than 30 */
/* SlidingWindow: performs a sliding window trimming approach: clips the read once the average quality within the window falls below a threshold: here for 5 bases, it is required an min average of 20 quality score */
/* MINLEN: Drop the read if it is below a specified length: here if it's below 50bp. Because too small reads are difficult to aligned */

"""
}


/*Test command Fastqc on the Trimmed fastq*/
/* We need to check if the quality of data will not decrease after the trimming, and we can continue to use it in other steps */

process QualityCTrimmed {
        module "FastQC"
        publishDir params.outdir, mode: 'copy'
	executor "lsf"

        input:
        set pair_id, file(treads) from TrimmomaticPaired

        output:
        set "*.zip", "*.html" into qualitycontrolRT
        set pair_id, file(treads) into TrimmomaticPaired2
"""
zcat ${treads}
fastqc ${treads}

"""
} 

/* End of the Trimming Steps: we have make all the necessary trimming, sorting and quality control to have better data for the alignment steps */

/* ALIGNMENT */
/* In the next steps we will aligned the reads data on a reference genome. For that we need to index the genome reference if the index not already exist. After that we can make alignment step with STAR and return bam files. We could plot the bam stats for check the quality of the bam files*/ 

/* Test command genome Index using Star: need reference genome and annotations files */

process buildIndex {
    executor "lsf"
    module "STAR-2.5.3a"
    publishDir params.outdir, mode: 'copy'

    input:
    file(genome) from genome_file   /* Input files: we need the reference genome in fasta format and the annotation file */
    file(annot) from annotation_file

    output:
    file('*.fa') into genome_index  /* Ouput file is the index reference genome in the genome_index channel */
    file(annot) into annotation_file2 /* We need to use annotation on next step, so as we have seen before we need to create another channel as we can't use twice the same channel */

    """
bash  STAR -genomeFastaFiles ${genome} -runThreadN 32 -limitGenomeGenerateRAM 120000000000 -sjdbOverhang 99 -sjdbGTFfile ${annot} -sjdbGTFtagExonParentTranscript Parent 
    """
}


/* Alignment with Star: Single Pass Mode: input are genome, genome index and fastqtrimmed data, output are BAM file*/

process algnment {
    executor "lsf"
    module "STAR"
    module "picard-tools"
    publishDir params.outdir, mode: 'copy'

    input:
    file index from genome_index
    set pair_id, file(treads) from TrimmomaticPaired2 /* We take the trimmed data as they are the data with the best quality score */

    output:
    set pair_id, file("*.bam") into bam /* We get back the bam files created by STAR alignment in the "bam" channel */
    
    """
    STAR -genomeDir ${index} -readFilesIn ${pair_id}.trimmomatic_1P.fq.gz ${pair_id}.trimmomatic_2P.fq.gz -chimSegmentMin 15 -chimJunctionOverhangMin 15 -alignMatesGapMax 20000 -alignIntronMax 20000 -outSAMtype BAM Unsorted -outSAMprimaryFlag AllBestScore -outSAMstrandField intronMotif -outFilterIntronMotifs RemoveNoncanonicalUnannotated -quantMode GeneCounts -outStd BAM_Unsorted -outQSconversionAdd -31
 """
}


/*Plot the BAM stats*/

process plotStat {
	executor "lsf"
	module "samtools"
	publishDir params.outdir, mode: 'copy'

	input:
	set pair_id, file(bamf) from bam
	
	output:
	file("*.stats") into statistics
	file("*.graph") into graphics
	set pair_idm file(bamf) into bam2
"""
samtools stats ${bam} > ${pair_id}.stats  /*Give some stats about the bam files*/
plot-bamstats ${pair_id}.stats > ${pair_id}.graph /*Give some graphics about the bam-stats*/
"""
}

/* End of the Alignment Steps */
	

/*Test Differential Expression Workflow*/
/* The next steps allow to see if there are differential expression between the expressed genes with some statistics */

/* HTSeq Count: get counts for differential expression */

process DiffCount {
	executor "lsf"
	publishDir params.outdir, mode: 'copy'
	module "HTSeq/0.7.2"

	input:
	file annot from annotation_file2
	set pair_id, file(bamf) from bam2
	
	output:
	file("*_HTseq.counts") into CountsDiff
"""
htseq-count -f bam -r pos -s no -i Parent ${bamf} ${annot} > ${pair_id}_HTseq.counts
"""
}

/*DESeq2: determine the statistically differentially expressed genes */
/*
process StatDiffExpr {
	library('ggplot2')
	library("gplots")
	library("vsn")
	library('RColorBrewer')
	library("DESeq2")
}
*/
	

/* Test Assembles the transcript by using the "cufflinks" tool */
/*
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'
    module "cufflinks"

    input:
    set pair_id, file(reads) from PairChannel(params.reads)
    file(bam_file) from bam

    output:
    file('transcript_*.gtf') into transcripts

    """
    cufflinks --no-update-check -q -p $task.cpus ${bam_file}
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}
*/

 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
