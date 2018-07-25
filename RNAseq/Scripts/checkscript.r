#!/usr/bin/env Rscript


rm(list=ls())                                        # remove all the objects from the R session
library(optparse)
library(crayon)

loadTargetFile <- function(targetFile, varInt, condRef, batch){
  target <- read.table(targetFile, header=TRUE, sep="\t", na.strings="")
  if (!I(varInt %in% names(target))) stop(paste(red("The factor of interest", varInt, "is not in the target file")))
  if (!is.null(batch) && !I(batch %in% names(target))) stop(paste(red("The batch effect", batch, "is not in the target file"))) 
  target[,varInt] <- as.factor(target[,varInt])
  if (!I(condRef %in% as.character(target[,varInt]))) stop(paste(red("The reference level", condRef, "is not a level of the factor of interest")))
  target[,varInt] <- relevel(target[,varInt],ref=condRef)
  target <- target[order(target[,varInt]),]
  rownames(target) <- as.character(target[,1])
  # check if varInt contains replicates
  if (min(table(target[,varInt]))<2) stop(paste("The factor of interest", varInt, "has a level without replicates"))
  # check if NA in the target
  if (any(is.na(cbind(target[,c(varInt, batch)], target[,1:2])))) stop("NA are present in the target file")
  # warning message if batch is numeric
  if (!is.null(batch) && is.numeric(target[,batch])) warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  if (any(grepl("[[:punct:]]", as.character(target[,varInt])))) stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
  cat("Target file:\n")
  print(target)
  return(target)
}

option_list <- list( make_option(c("-t", "--targetFile"), default="target.txt",dest="targetFile",help="path to the design/target file [default: %default]."), make_option(c("-v", "--varInt"), default="group", dest="varInt", help="factor of interest [default: %default]"),make_option(c("-b", "--batch"),default=NULL,dest="batch",help="blocking factor [default: %default] or \"batch\" for example"),make_option(c("-c", "--condRef"),default="WT",dest="condRef",help="reference biological condition [default: %default]"))

parser <- OptionParser(usage="usage: %prog [options]",
                                           option_list=option_list,
                                           description="Compare two or more biological conditions in a RNA-Seq framework with DESeq2.",
                                           epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

workDir <- getwd()
targetFile <- opt$targetFile 
varInt <- opt$varInt                               
condRef <- opt$condRef                               
batch <- opt$batch 

target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

