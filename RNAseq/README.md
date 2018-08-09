

# RNA-seq pipeline

[![nextflow]( https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)

Author: Carole Smerilli

## Overall

This pipeline was developed based on Nextflow language, inspired on RNA-Seq Analysis Differential Expression Notebook [Jupyter]( https://jupyterhub.powerplant.pfr.co.nz/user/cfpcxs/notebooks/cfpcxs/RNASeq_community/best_practice_drafts/RNA_Seq_Analysis_DE_Best_Practices_Rev_November_2017.ipynb) written by Charles David. This pipeline is able to apply the main steps of **RNA-seq analysis**, like **Quality Control**, **Trimming**, **Alignment**, **Counts** and **Differential Expression** on fastQ files which contains **Paired-end reads**, with a **reference genome**. It contains more steps of Quality Control to check the quality of data during all the trimming steps, is able to edit a database which contains counts, and take in account potential human error. The nextflow.config file allows the users to change a lot of parameter and choose options or tools they prefer to analyse data (for example DESeq and/or EdgeR…). 

## Overview of the structure





![image](https://user-images.githubusercontent.com/39930402/43873893-be9e1634-9bdc-11e8-87c0-0e8974690c18.png)









  ![flowchart](https://user-images.githubusercontent.com/39930402/43873830-75fb91f4-9bdc-11e8-9ce9-f78635de5911.png)
















  



  












## Table of Contents
- [Mandatory input files](#mandatory-input-files)
- [Installation](#installation)
- [Configuration](#configuration-file)
- [Conditions file](#conditions-file)
- [Database](#database)
- [Error/Crash](#error-job-killed-crash)

## Mandatory input files

There are some necessary data that pipeline need to know. Please check that you have these files before try to use it.

-	Paired-End read on **FastQ files** format issued of Illumina sequencing (If you have several files due to different lanes you need to merge the files before use this pipeline)
-	**Genome reference**
-	**Annotation** (GTF/GFF format) to build the index of genome reference
-	**Conditions file** in a specific format (see the documentation part below)
-	Genome contamination index, database (they are already provided in the default configuration file but you could change if you want)

## Installation
This pipeline is implemented with Nextflow language, so it requires Nextflow to be installed. Inside the pipeline the use of different language and tools requires Bash, Python, R, SQL to be installed. 
Clone *main.nf*, *nextflow.config* and the *Scripts* directory in the same environment than Nextflow. (Command `git clone https://...`)

## Configuration file 

### Parameters

The *nextflow.config* file contains all the variable you could change, mandatory variable and optional variables. You need to open it (command `nano nextflow.config` ) to change variable.
See below an overview of the configuration file. All the lines after "=" sign should/could be replace by your parameters, as you can see there is an explanation for each parameter in comment. 



	 standard{
	

		params{
			email='carole.smerilli@plantandfood.co.nz'
	

			namefastqfileR1= '_R1.fastq'    /* You need to define the structure of the name of fastq file: for example DA8_sub_R1.fastq, write '_R1.fastq', if it is DQ8_R1.fq, write '_R1.fq'. So write without the ID */ 
			namefastqfileR2= '_R2.fastq'    /* Paired read so for Read and Forward read */
			species= 'Potatoes'            /* Please inform the species of your data */
			conditionsfile= "${baseDir}/targets.txt"  /* Please inform the location of the file which contains a table with the name of the file, ID and conditions. */
			reads= "/workspace/cfpcxs/PFRAutomatedWorkflows/RNAseq/Potatoes/*_R{1,2}.fastq"  /* Inform the location of your fastq files and put the right structure '_R{1,2}' to collect them. */
			annot= "/workspace/cfpcxs/PFRAutomatedWorkflows/RNAseq/Potatoes/PGSC_DM_V403_fixed_representative_genes.gff" /* Defines the path for the annotation file */
			genome= "/workspace/cfpcxs/PFRAutomatedWorkflows/RNAseq/Potatoes/PGSC_DM_v4.03_pseudomolecules_ALL.fasta"  /* Defines the part for the reference genome file (fasta format) */
			     }
       }
			

### Profiles
It contains different profiles which are all explains in the file and give you choice to use an already existing Index, sample of your data, DESeq and/or EdgeR. To run this pipeline you need to inform profiles in a cumulative way. 
For example: 

`./nextflow run main.nf -profile standard,IndexE,sampleN,EdgeR `

This command run the pipeline with the profile standard, IndexE, sampleN and EdgeR :
Be careful to always inform at least one of each profiles, that is to say IndexN if you don't want IndeE, sample if you don't want sampleN.... 

**Standard profile:** contains all the basic parameter like the location of the Fastq files, conditions file, species, the options of each tools (trimmomatic, STAR…), the output directory.... You always need to inform this profile in the command line!

**IndexN profile:**  if you have not the index of the reference genome for these data you need to build it, so you need to use the IndexN profile.

**IndexE profile:** if you already have the index, put the address of his location in the IndexE profile, and call IndexE when you run the pipeline, it will save time.

**sample profile:** if you want to test the pipeline only on sample of your data, more small than your data, call the sample profile when you run the pipeline.

**sampleN profile:** if you want to run the pipeline on your full/complete data, call the sampleN profile when you run the pipeline.

**DESeq profile:** if you want to make the Differential Expression statistics with DESeq packages, you could change parameter in profile Standard for common options and in DESeq profile for specific option. You could cumulate DESeq with EdgeR.

**EdgeR profile:** if you want to make the Differential Expression statistics with EdgeR packages, you could change parameter in profile Standard for common options and in EdgeR profile for specific option. You could cumulate EdgeR with DESeq.

## Conditions file

You need to inform the targets file, which will be used in the Differential Expression analysis (DESeq, EdgeR). For this pipeline it is necessary to respect a specific format. Please follow these indications.
 
It must be composed at least of three columns: the first one is the ID column (contain the ID of sample), the second is the file column (contain the name of the file, please note that for the moment it is necessary to replace “R*.fastq” by “HTseq.counts” but it could be change thereafter) and the third one is the conditions column (contain the different conditions of the experiment). Each columns need to be separated by a tabular space. 

Depending of the experiment you could have a *batch effect* for example, an additional condition depending of the day or of tissues. It is necessary that you add a column to the target file which give information about this addable conditions, you could call this column the name you want (*day* for example).

Please take in account that you will need to check in the configuration file of nextflow (*nextflow.config*), in the **Standard profile**,  part *Parameters for the Differential Expression statistics*, if the variables correspond to the name of the columns on your targets file.  

For example: 

```
id	file	group	rep
8HPT	8HPT_sub_HTseq.counts	hot psyllid	rep1
6HPT	6HPT_sub_HTseq.counts	hot psyllid	rep2
11HPT	11HPT_sub_HTseq.counts	hot psyllid	rep3
2CT	2CT_sub_HTseq.counts	control	rep1
4CT	4CT_sub_HTseq.counts	control	rep2
12CT	12CT_sub_HTseq.counts	control	rep3
```

Here the condition are ‘Control’ and ‘Hot psyllid’. The name of the conditions column is “group” and there is not batch effect. So in the configuration file (*nextflow.config*) you need to inform the variables like that.

```
factorOfInterest="--varInt 'group'"    
ReferenceCondition="--condRef 'control'"  
BlockingFactor="--batch NULL”
```
If you have a batch effect, for example a day condition, you need to create another column, for example:

```
id	file	group	day	rep
8HPT	8HPT_sub_HTseq.counts	hot psyllid	0	rep1
6HPT	6HPT_sub_HTseq.counts	hot psyllid	3	rep2
11HPT	11HPT_sub_HTseq.counts	hot psyllid	7	rep3
2CT	2CT_sub_HTseq.counts	control	0	rep1	
4CT	4CT_sub_HTseq.counts	control	3	rep2
12CT	12CT_sub_HTseq.counts	control	7	rep3
```

In this example the experiment is running during different period of time (0,3 or 7 days). So here you have another column which inform these different period, the name of this column is *day* for example.
If you want to take in account the time effect you need to inform the batch effect in the configuration file, in giving the name of the column:
```
factorOfInterest="--varInt 'group'"    
ReferenceCondition="--condRef 'control'"  
BlockingFactor="--batch ‘day’"
```
But if you don’t want to take the batch effect in account you could just let the batch effect *NULL*:

```
BlockingFactor="--batch NULL”
```
This pipeline allow to run a general differential analysis and give you a first global approach of statistics analysis. If you want more specific analysis you probably need to make theses analysis by yourself for the moment. 

## Database
This pipeline edit a database (see below) that allow every people on Plant and Food to check and retrieve information about each genes and the number of counts for each sample in easy way by using database management. This database contains 6 columns: ID, Experiment’s ID, ID of sample, name of gene, number of counts and species. 

![image](https://user-images.githubusercontent.com/39930402/43873801-50e9c976-9bdc-11e8-991e-3e1f4b41ce89.png)

The original address of this database is: *cfphxd_rnaseq_counts*. To access to the database you need the username: *cfphxd_rnaseq*. Password: *yIPGRMPrDjtgUTeu*. 
The name of the table is “databaseT”.
 
If you want to change this you need to create another table. You need to make a connection with dbaas, in the shell, write this kind of command (replace the name of user, password and name of the database):

`mysql -h database.powerplant.pfr.co.nz -p password  -u username -D databasename`

Once the connection to dbaas is working you could write this kind of command in the SQL environment (you could replace the name of the table ‘databaseT’):

`CREATE TABLE databaseT ( ID int NOT NULL AUTO_INCREMENT, ExperimentID varchar(20), SampleID varchar(25) NOT NULL, Gene varchar(30) NOT NULL, Counts INT, Species varchar(20) NOT NULL, primary key (ID) );`

## Error, Job killed, Crash

### Resume option:

 If the pipeline has crashed because of a lack of memory or too busy servor, you could use the **resume option**, this option allow to start at the steps where the pipeline crashed and save the steps/processes already run by the pipeline before. So it’s avoid to lose time to rerun all the previous steps. You could use this option like that: 
 
`./nextflow run -resume main.nf profiles standard,IndexE,sampleN,EdgeR`

### Error:

If an error occured, the pipeline is going to stop and return an error message which explain the reason. Some process have the purpose to return an error message like the quality control test if the crucial step fail the test, or also the process which check the conditions file if its in a bad format. Please pay attention of the error message you get.





