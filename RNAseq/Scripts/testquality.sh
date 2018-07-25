#!/bin/bash

# This scrip is used in all QualityTest Process. It checks if the second line of the summary file of FastQC report Pass or Fail: check if the Per Base Quality is ok or not. If Fail, all the pipeline script is stopped and an error message is return to ask the user to check its data and remove the fail data if he/she can.


# $1 first parameter: zip file of fastqc report
# $2 second parameter: pair_id.nomfichier

unzip $1

function traitement() {

var=$(sed -n "2p" $1/summary.txt | grep "FAIL" | awk {'print $1'})

if [ -z $var ];
	then value="1"
        else value="2"
fi
return ${value}
}

traitement $2
echo $value

if [ $value = "2" ];
	then echo "The fastq file $2 has too bad Per Base Quality, please check your fastq file and remove it if you can"; exit 1
fi
