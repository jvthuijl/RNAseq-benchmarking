# Simulate reads assembled transcript
# usage: Rscript base_expression.r fileIn fileOut N

# N = reads per transcript
# source("http://bioconductor.org/biocLite.R")
# biocLite("polyester")
# biocLite("Rsubread")
library(polyester)
library(Biostrings)

# reset variables
rm(list=ls(all=TRUE))

args<-commandArgs(TRUE)
set.seed(42)

fastapath = args[1]
fasta = readDNAStringSet(fastapath)

len = length(fasta)
print(len)

# base of 30
fold_changes = matrix(30, nrow=len)

paste("Simulating", args[3], " reads.")
simulate_experiment(fastapath, paired = TRUE, fold_changes = fold_changes,
                    reads_per_transcript = as.integer(args[3]), num_reps = c(20), 
                    outdir=paste(args[2], 'base', args[3], sep = "_"), seed = 42 )
warnings()
