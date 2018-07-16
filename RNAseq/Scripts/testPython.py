#!/usr/bin/python


# This Python script is used in the CountsTable process. It take the HTSeq_Counts files in input data and return a unique file in txt format with a table :  the headers (horizontal) of this table are the sample ID of the files, and the first column are the name of the genes, so you could see for each gene, how many counts each sample have 




import glob

ColumnHeader=[] # This list will contains all the header columns like the sample ID

genes={} #Dictionnary for each gene like genes={ gene1= [0,1,0,3], gene2=[2,0,1,1], ...}

outfile= open('merged_counts.txt','w') #Create the output file with the table

for file in glob.glob("*_HTseq.C.counts"):  #For file of counts data in the actual directory, the temporary directory where the process is running
	ColumnHeader= ColumnHeader + [file.split('_')[0]]  #Take back the Sample_ID of each files, and insert in the list of Column Header
	counts_f= open(file)  
	
	for line in counts_f: #For each line of this file, so for each gene and its counts
		gene,counts= line.strip().split('\t') #Take back the gene and its associate counts
		if gene not in genes:  #If the gene doesn't exist in the dictionnary
			genes[gene]=[counts] #Create the appropriate key with counts associate
		else:
			genes[gene].append(counts) #If the gene already exist in the dictionnary, just complete the list of counts
counts_f.close()

print ColumnHeader

outfile.write('Gene\t'+'\t'.join(ColumnHeader)+'\n') #Write the first lines, with headers and names of columns

for gene in genes:   #Run through the genes dictionnary and for each gene
	count= genes[gene] #Get back the counts list associated with that gene
	out_counts= gene+'\t'+'\t'.join(count)+'\n' #Formatting data for the table
	outfile.write(out_counts) #Write in the table

outfile.close()

print "Merged table is 'merged_counts.txt'"
		
		
		


