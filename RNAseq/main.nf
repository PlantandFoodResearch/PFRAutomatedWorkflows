
/*
 *   RNA -seq Pipeline 
 */



/* REQUIREMENTS */
/* To run this  pipeline you need: Fastq files Paired-end Reads, reference genome in fasta format, annotation in gff format, a database contamination genome in fasta format and their index in db format. You need to define an output directory in which you will find all the output data (params.outdir), you need to know the structure of the name of the fastq files like ID_sub_R1.fastq or ID.R1.fq or ID.F1.fastq. You need to copy the sh script in your environment, it will be useful for the test of quality */


params.outdir = 'results54'   /* Defines an output directory where you will be able to see the output files/data of this pipeline */
 


if(params.A=='0'){

        Channel.fromFilePairs(params.reads).set{ReadChannelC}}

else {
        Channel.empty().set{ReadChannelC}
}

Channel.fromFilePairs(params.reads).ifEmpty { error "Cannot find any reads matching: ${params.reads}" }.set{ReadChannel}

 

/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)



/* process check allow to check if the conditions file is well formated. As it's the most sensitive and risky part of the pipeline and subject to human mistake, it's more useful to check at the beginning of the pipeline rather than return an error at the end of all the pipeline. */

process check {

	publishDir params.outdir
	conda "/workspace/appscratch/miniconda/cfpcxs_sartools"
	
"""
${params.scriptCheck} --targetFile ${params.conditionsfile} ${params.factorOfInterest} ${params.ReferenceCondition} ${params.BlockingFactor} 
"""

}



/* First step */

/*Test on small sample of real data: test if the pipeline works*/ 

process sample {
       module "seqtk"
       publishDir params.outdir, mode:'copy', overwrite:'true'
       

        input:
        set pair_id, file(reads) from ReadChannel

        output:
        set pair_id, file("*.fq") into Sample
	
	when:
	params.A=='1'
"""
seqtk sample -s100 ${pair_id}${params.namefastqfileR1} ${params.size} > ${pair_id}.R1.fq
seqtk sample -s100 ${pair_id}${params.namefastqfileR2} ${params.size} > ${pair_id}.R2.fq
"""
}


/* Test command Fastqc: take Fastq files which contains the reads and return Quality control report of these reads in html and zip format. There is the first check of the quality of data juste after the sequencing. */

process qualitycontrol {   /* Copy the return data in the output directory: allow to see the data files  */
	label 'QC' 
	publishDir params.outdir1, mode: 'copy', overwrite:'true'

        input:
        set pair_id, file(reads) from Sample.mix(ReadChannelC) /*Sample contains small samples of the reads data: the pair_id is an ID in common for each forward and reverse read, for example if you have Kiwitestdata1.R1.fq, Kiwitestdata1.R2.fq, the pair_id will be "Kiwitestdata1", file(reads) just take the corresponding reads of the Sample channel */

        output:
        set "*.html", "*.zip" into qualitycontrolR /* Get back the QCreports returned by the fastqc command in the channel QualitycontrolR */
        set pair_id, file(reads) into Sample2  /* In Nextflow, you can't use the same channel twice, so as we need to use the reads containing in the channel sample in other process, we need to stock the reads data in other channel, for example: Sample2 */

/* fastqc is a command which take the reads and return QCreports*/

"""
fastqc ${reads}  
"""
}

/* Test if the quality of the reads is good enough to pass at the next step: use the "testscript.sh": check if the quality per base is ok. If not: return an error with a message tell the user to check the fastq fil which FAIL the quality test and remove if possible. */

process Testqualitycontrol {
         publishDir params.outdir, mode: 'copy'
         

        input:
        file(readzip) from qualitycontrolR
        set pair_id, file(reads) from Sample2

        output:
        file(readzip) into qualitycontrolR2
        set pair_id, file(reads) into Sample3



"""
bash ${params.scriptTQC} ${pair_id}.R1_fastqc.zip ${pair_id}.R1_fastqc
bash ${params.scriptTQC} ${pair_id}.R2_fastqc.zip ${pair_id}.R2_fastqc

"""

}


/* TRIMMING steps */
/* In the next steps, we will check the quality of reads and if necessary make some trimming to enhance the quality, before to go on the alignment step */

/* Test command SortMeRNA: it allows to remove the contaminated reads. We need to merge the reads before using the SortMeRNA filter, with the bash script "merge-paired-reads" we obtain merged/interleaved fastq files. With these files we make the SortMeRNA and getabck only the data which are not contaminated. We need to unmerge the data with the script "unmerge-paired-reads" to obtain fastq files unmerged in the aim to use them in next steps*/

process SortMeRNA {
	publishDir params.outdir, mode: 'copy', overwrite:'true'

	module "sortmerna/2.1"
	

	input:
	set pair_id, file(reads) from Sample3

	output:
	file('*.MERGED.fastq') into MergedData
	file('*_filtered.fastq') into FINTERLEAVEDFiltered
	set pair_id, file('*_R*.fq') into FQFiltered  /* Here, we get back the final unmerged data sorted in the FQFiltered channel: reads file in fastq format without the contaminated data */


/* Probably need to gunzip with zcat <(zcat filename.gz) */
/* Aligned the reads data with different contamination genome inquired in the "--ref", the variable  "{SORTMERNADB}" need to be inquired by the user on the top of the pipeline by the user. The  aligned data are the contaminated data so we remove them and just get back the data which were aligned data are the contaminated data so we remove them and just get back the data which were  not aligned with the contamination genome so the data ouput in " --other" */


"""
bash merge-paired-reads.sh ${pair_id}${params.subd1} ${pair_id}${params.subd2} ${pair_id}.MERGED.fastq 
sortmerna --ref ${params.SORTMERNADB} --reads ${pair_id}.MERGED.fastq --paired_in --fastx --aligned ${pair_id}_rRNA --other ${pair_id}_filtered    
bash unmerge-paired-reads.sh ${pair_id}_filtered.fastq ${pair_id}f_R1.fq ${pair_id}f_R2.fq  
"""
}



/*Test command Fastqc on the filtered fastq*/
/* We need to check if the quality of data will not decrease after the sorting, and we can continue to use it in other steps */

process QualityCFiltered {
        label 'QC'
	publishDir params.outdir1, mode: 'copy', overwrite:'true'


        input:
        set pair_id, file(freads) from FQFiltered

        output:
        set "*.zip", "*.html" into qualitycontrolRF
        set pair_id, file(freads) into FQFiltered2
"""
fastqc ${freads}
"""
} 

process Testqualitycontrol2 {
         publishDir params.outdir, mode: 'copy'
         

        input:
        file(readzip) from qualitycontrolRF
        set pair_id, file(reads) from FQFiltered2

        output:
        file(readzip) into qualitycontrolRF2
        set pair_id, file(reads) into FQFiltered3



"""
bash ${params.scriptTQC} ${pair_id}f_R1_fastqc.zip ${pair_id}f_R1_fastqc
bash ${params.scriptTQC} ${pair_id}f_R2_fastqc.zip ${pair_id}f_R2_fastqc

"""

}


/*Test command Trimmomatic*/
/* Sometimes we need to make some trimming in the data for inhance the quality of data, here we use the trimmomatic script */

process Trimmomatic {
        module "java"
        module "Trimmomatic-0.36"
        publishDir params.outdir2, mode: 'copy', overwrite:'true'

	

        input:
        set pair_id, file(fqreads) from FQFiltered3

        output:
        set pair_id, file('*trimmomatic_*P.fq.gz') into TrimmomaticPaired /* We get back the data which are correctly trimmed (have survived of the trimming): both of the reads (forward and reverse) have survived of the trimming, ine the TrimmomaticPaired channel */
        set pair_id, file('*trimmomatic_*U.fq.gz') into TrimmomaticUnpaired /*  We get back the data which are not correctly trimmed: at least one of the paired-reads (forward or reverse) or the two have be removed by the trimming, ine the TrimmomaticUnPaired channel */

"""
java -jar /software/bioinformatics/Trimmomatic-0.36/trimmomatic.jar PE ${pair_id}f_R1.fq  ${pair_id}f_R2.fq -baseout ${pair_id}.trimmomatic.fq.gz ILLUMINACLIP:${params.Adapterfile}${params.Pclip} SLIDINGWINDOW${params.SlidWindow} MINLEN${params.MinL}
"""
/* We need to indicate if it's pair-ended reads with "PE", just after that we put the reads data, in the -baseout it will return the ouput file for each pair ended read so the Paired and UnPaired Trimmomatic for forward and reverse reads so 2 file for Paired (R1Paired, R2Paired) and 2 files for UnPaired (R1UnPaired, R2UnPaired) */
/* Illuminaclip: find and remove Illumina adapters which could contaminated the reads data during the sequencing step: it's need a fasta file containing illumina adapters (here we inquire the pathway of this file), the "30" is for clipping when the score is lower than 30 */
/* SlidingWindow: performs a sliding window trimming approach: clips the read once the average quality within the window falls below a threshold: here for 5 bases, it is required an min average of 20 quality score */
/* MINLEN: Drop the read if it is below a specified length: here if it's below 50bp. Because too small reads are difficult to aligned */


}


/*Test command Fastqc on the Trimmed fastq*/
/* We need to check if the quality of data will not decrease after the trimming, and we can continue to use it in other steps */

process QualityCTrimmed {
        label 'QC'
	publishDir params.outdir1, mode: 'copy', overwrite:'true'


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
process Testqualitycontrol3 {
         publishDir params.outdir, mode: 'copy'
         

        input:
        file(readzip) from qualitycontrolRT
        set pair_id, file(reads) from TrimmomaticPaired2

        output:
        file(readzip) into qualitycontrolRT2
        set pair_id, file(reads) into TrimmomaticPaired3



"""
bash ${params.scriptTQC} ${pair_id}.trimmomatic_1P_fastqc.zip ${pair_id}.trimmomatic_1P_fastqc
bash ${params.scriptTQC} ${pair_id}.trimmomatic_2P_fastqc.zip ${pair_id}.trimmomatic_2P_fastqc

"""

}
/* End of the Trimming Steps: we have make all the necessary trimming, sorting and quality control to have better data for the alignment steps */




/* ALIGNMENT steps */
/* In the next steps we will aligned the reads data on a reference genome. For that we need to index the genome reference if the index not already exist. After that we can make alignment step with STAR and return bam files. We could plot the bam stats for check the quality of the bam files*/ 



/* Test command genome Index using Star: need reference genome and annotations files */
/* process Index check if the genome's index already exist and use it, if not it build the index. This conditional process avoid to lose time to build index ifit's  already existing. */

process Index {
    
    module "STAR/2.5.3a"
    publishDir params.outdir5, mode: 'copy', overwrite:'true'


    input:
    file(genome) from genome_file   /* Input files: we need the reference genome in fasta format and the annotation file */
    file(annot) from annotation_file

    output:
    file("indexDirB") into indexC  /* Ouput file is the index reference genome in the genome_index channel */
    file(annot) into annotation_file2 /* We need to use annotation on next step, so as we have seen before we need to create another channel as we can't use twice the same channel */
    

script:

if( params.skip=='1')
"""
mkdir indexDirB

STAR --runMode genomeGenerate --genomeDir indexDirB --genomeFastaFiles ${genome} --sjdbGTFfile ${annot} ${params.StarAddition} 
"""

else
"""

ln -s ${params.index} indexDirB
"""
 } 

/* Alignment with Star: Single Pass Mode: input are genome, genome index and fastqtrimmed data, output are BAM file*/

process alignment {
    module "STAR/2.5.3a"
    module "picard-tools"
    publishDir params.outdir3, mode: 'copy', overwrite:'true'


    input:

    file(index1) from indexC
    set pair_id, file(treads) from TrimmomaticPaired3 /* We take the trimmed data as they are the data with the best quality score */

    output:
    set pair_id, file("BAMD.*") into bam /* We get back the bam files created by STAR alignment in the "bam" channel */
    
"""
mkdir BAMD.${pair_id}
STAR --genomeDir ${index1} --outFileNamePrefix BAMD.${pair_id}/${pair_id}_ --readFilesCommand zcat --readFilesIn ${pair_id}.trimmomatic_1P.fq.gz ${pair_id}.trimmomatic_2P.fq.gz --chimSegmentMin ${params.SegmMin} --chimJunctionOverhangMin ${params.JunctMin} --alignMatesGapMax ${params.GapMax} --alignIntronMax ${params.IntronMax} --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode GeneCounts --outQSconversionAdd ${params.outQSconvAdd}
"""
}


/*Plot the BAM stats*/

process plotStat {
	module "samtools"
	publishDir params.outdir4, mode: 'copy', overwrite:'true'


	input:
	set pair_id, file(bamf) from bam
	
	output:
	file("*.stats") into statistics
	file("statisticsG.*") into graphic
	set pair_id, file(bamf) into bam2
"""
mkdir statisticsG.${pair_id}
samtools stats ${bamf}/${pair_id}_Aligned.out.bam > ${pair_id}.stats  
plot-bamstats ${pair_id}.stats -p statisticsG.${pair_id}/${pair_id}
"""
/*Give some stats about the bam files*/
/*Give some graphics about the bam-stats*/

}

/* End of the Alignment Steps */
 	


/*DIFFERENTIAL EXPRESSION Workflow steps*/

/* The next steps allow to see if there are differential expression between the expressed genes with some statistics */

/* HTSeq Count: get counts for differential expression and database creation */
/* Return HTseqcounts files, remove the last five lines of the HTSeqcounts (only for insert data, the output Counts files are full) and insert data of this alignment in a database: use the counts, gene, pair_id of the HTseq files, the species need to be inform by user. The loop while use these data in paramater of the bash script "testSQL3 which make SQL command ot insert data in a database. The database need to be provided / created by user as well as the user name, code, servor  and database name. */

process DiffCount {
	publishDir params.outdir6, mode: 'copy', overwrite:'true'

	module "HTSeq/0.7.2"
	module "mysql.connector"

	input:
	file annot from annotation_file2
	set pair_id, file(bamf) from bam2
	
	output:
	set pair_id, file("*_HTseq.counts") into CountsDiff
	set pair_id, file("*_HTseq.C.counts") into CountsCDiff

"""
htseq-count -f bam -r pos -s no -i Parent ${bamf}/${pair_id}_Aligned.out.bam ${annot} > ${pair_id}_HTseq.counts
head -n -5 ${pair_id}_HTseq.counts > ${pair_id}_HTseq.C.counts
while read LINE; do bash ${params.scriptSQL} \$LINE ${params.species} ${pair_id} ${params.IDexp} >> buffer.sql; done < ${pair_id}_HTseq.C.counts
mysql -hdatabase.powerplant.pfr.co.nz -p${params.passwordDB} -u${params.userDB} -D${params.nameDB} < buffer.sql 

"""
}
/* SQL Server buffer is a place in system memory that is used for caching table and indexx data pages as they are modified or read from disk */

/* process TableCounts created a text file with table containing all the counts of each sample for each genes, the same gene is not in replicats so it's more readable. */
/* So it's using the testPython script: a Python language script to build this table from theHTseq.C.counts data. The use of ".collect()" allow to collect all the channel of HTseqcounts in an unique channel, so the process run only one time, and not 6 times. */

process TableCounts {
	
	publishDir params.outdir6, mode: 'copy', overwrite:'true'


	input:
	file(CountsCF) from CountsCDiff.collect() 

	output:
	file("merged_counts.txt") into MergedC
	
"""
ls -l
python ${params.scriptTabC} 
"""
}

/* process StatDE is  the statistics analysis step for differential expression workflow.  It is using SARTools, a script which allow the use of DESEq2 and/or EdgeR packages and return all graphics, tables and necessary files to measure differential expression of transcripts in different conditions.*/
/* This process need the HTSeq counts files, the conditions file and some optional variables that user could change in nextflow.config. */

process StatDE {
        conda "/workspace/appscratch/miniconda/cfpcxs_sartools"
        publishDir params.outdir7, mode:'copy', overwrite:'true'


	input:
	file(HTSeq) from CountsDiff.collect()
        
	output:
        file("DESeq_Output") optional true into DESEQreport
	file("EdgeR_Output") optional true into EdgeRreport
	file("figures") into Figure
	file("tables") into Table

script:

if(params.DESeq=='1' & params.EdgeR=='1')
"""
${params.scriptDESeq} --targetFile ${params.conditionsfile} --rawDir . ${params.featuresToRemove} ${params.factorOfInterest} ${params.ReferenceCondition} ${params.BlockingFactor} > DESeq_Output
${params.scriptEdgeR} --targetFile ${params.conditionsfile} --rawDir . ${params.featuresToRemove} ${params.factorOfInterest} ${params.ReferenceCondition} ${params.BlockingFactor} > EdgeR_Output
"""

else if(params.EdgeR=='1')
"""
${params.scriptEdgeR} --targetFile ${params.conditionsfile} --rawDir . ${params.featuresToRemove} ${params.factorOfInterest} ${params.ReferenceCondition} ${params.BlockingFactor} > EdgeR_Output
"""

else if(params.DESeq=='1')
"""
${params.scriptDESeq} --targetFile ${params.conditionsfile} --rawDir . ${params.featuresToRemove} ${params.factorOfInterest} ${params.ReferenceCondition} ${params.BlockingFactor} > DESeq_Output
"""


}



/*Notification message: send an e-mail when the worklow execution is completed or an report error if the pipeline has been stopped. */

workflow.onComplete {
    def subject = 'My pipeline execution'
    def recipient = params.email

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
} 
