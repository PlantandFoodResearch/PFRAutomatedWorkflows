
#!/bin/bash

# This script allows to insert the data of the HTseqCOunts files in a database, with a ExpID automatically generated

# Below you could find the necessary parameter for this script:
# $1 = gene
# $2 = counts 
# $3 = species
# $4 = SampleID  = pair_id

gene=$1
count=$2
species=$3
sampleid=$4


echo "INSERT INTO databaseT (SampleID, Gene, Counts, Species) VALUES ( '$sampleid', '$gene' , '$count' , '$species') ; "





