# Simulate reads transcript differential expression
# usage: Rscript diff_expression.r fileIn fileOut N

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

stepsize = 2000

fastapath = args[1]
fasta = readDNAStringSet(fastapath)

numberOfTranscripts = length(fasta)
print(numberOfTranscripts)

# fixed levels of expression
levels = matrix()
levels[1:6] = 8
levels[6:11] = 6
levels[11:16] = 4
levels[16:21] = 2
levels[21:26] = 1.8
levels[26:31] = 1.6
levels[31:36] = 1.4
levels[36:41] = 1.2
levels[41:46] = 0.90
levels[46:51] = 0.8
levels[51:56] = 0.6
levels[56:61] = 0.4
levels[61:65] = 0.125

# base of 30
normal = matrix(1, nrow = numberOfTranscripts)
diff = normal
diff_positions = sample(1:numberOfTranscripts, 65)
print(diff_positions)
diff[diff_positions] = sample(levels, 65)

fold_changes = cbind(normal, diff)
print(fold_changes[diff_positions])
write.table(fold_changes, paste(args[2], "fold_changes", args[3], ".txt", sep = "_"))

print(paste("Simulating", args[3], "reads per transcript.")) 
simulate_experiment(fastapath, paired = TRUE, fold_changes = fold_changes,
                    reads_per_transcript = as.integer(args[3]), num_reps=c(10, 10),
                    outdir=paste(args[2],'diff', args[3], sep = "_"), seed = 42)
warnings()
