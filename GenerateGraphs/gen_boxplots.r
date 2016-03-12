library(ggplot2)
library(data.table)
library(splitstackshape)
library(corrplot)

rm(list=ls(all=TRUE))
wd <- "/home/umcg-jvanthuijl/umcg-jvanthuijl/kallisto"
setwd(wd)

# function to retrieve the abundance levels for all experiments within a directory dir.
calcAbundance <- function(dir, nrows = 20){
  outputAbundance = data.frame()
    prevWd <- getwd()
  setwd(dir)
  #for each samplefolder
  for (i in dir(dir)) {
    print(i) 
    #read the abundance levels
    tmpAbundance <- read.table(paste0(i, '/abundance.tsv'), nrows=nrows,
                               header = TRUE)
    # if target_id does not yet exist in the dataframe, add this column.
    if (!("target_id" %in% names(outputAbundance)))  {
      outputAbundance <- data.frame(tmpAbundance$target_id)
      colnames(outputAbundance) <- "target_id"
    }
    # for each sample add a column
    outputAbundance[[i]] <- tmpAbundance$tpm
  }
    setwd(prevWd)
  # create a new dataframe, with the target_id and the mean of all samples
  outputAbundance <- data.frame(row.names = outputAbundance[,1],
                                mean = rowMeans(outputAbundance[,-1]),
                                sd = apply(outputAbundance[,-1],1, sd),
                                sample = i)
  outputAbundance
}

calcDiffAbundance <- function(listSamples, nrows = 20){
  outputAbundance = data.frame()
  prevWd <- getwd()
  setwd("/home/umcg-jvanthuijl/umcg-jvanthuijl/kallisto")
  #for each samplefolder
  for (i in listSamples) {
    print(i) 
    #read the abundance levels
    tmpAbundance <- read.table(paste0(i, '/abundance.tsv'), nrows=nrows,
                               header = TRUE)
    # if target_id does not yet exist in the dataframe, add this column.
    if (!("target_id" %in% names(outputAbundance)))  {
      outputAbundance <- data.frame(tmpAbundance$target_id)
      colnames(outputAbundance) <- "target_id"
    }
    # for each sample add a column
    outputAbundance[[i]] <- tmpAbundance$tpm
  }
  setwd(prevWd)
  # create a new dataframe, with the target_id and the mean of all samples
  outputAbundance <- data.frame(target_id = outputAbundance[,1],
                                mean = rowMeans(outputAbundance[,-1]),
                                sd = apply(outputAbundance[,-1],1, sd),
                                sample = i)
  outputAbundance
}

# set the samplegroups. Folder stucture:
# results
#   |-cdna
#     |-10x
#     |  |- sample001
#     |  |- sample002 ...
#     |-20x
#        |- sample001
#        |- sample002 ...
samplegroups  <- data.frame()
abundances <- list()
abundances_diff <- list()
corrr <- list()

kmers=c("k11", "k21", "k31")
coverage=c("10x", "15x", "20x")

samplegroups = paste("cdna", rep(coverage, each=3), kmers, sep="/")
sampleprefixes = paste("cdna", rep(coverage, each=3), kmers, sep="_")
# pdf( "CoordPlots.pdf", width = 10, height = 10 )
# first = TRUE
# for (group in samplegroups) {
#   #read the abundance levels
#   abundance = calcAbundance(paste(wd, "base", group, sep = "/"), 400)
#   # if target_id does not yet exist in the dataframe, add this column.
#   if (first)  {
#     first = FALSE
#     print(head(abundance))
#     final <- data.frame(row.names = rownames(abundance))
#     #colnames(final) <- "target_id"
#   }
#   # for each sample add a column
#   final <- data.frame(final, new = abundance$mean)
#   names(final)[names(final) == 'new'] <- gsub("/", ".", group)
# }
  
# print(head(final))

setwd("/home/umcg-jvanthuijl/umcg-jvanthuijl/source/R/gen_graphs/")

# write.table(final, "illumina_probe.csv", col.names=TRUE, sep=",")

# set lower and upper threshold for the graph.
threshold = c(0.4, 1.3)

first = TRUE
for (group in samplegroups) {
  #read the abundance levels
  abundance = calcAbundance(paste(wd, "base", group, sep = "/"), 20000)
  # if target_id does not yet exist in the dataframe, add this column.
  if (first)  {
    first = FALSE
    print(head(abundance))
    means <- data.frame(target_id = rownames(abundance))
  #colnames(final) <- "target_id"
  }
  # for each sample add a column
  means <- data.frame(means, new = abundance$mean)
  names(means)[names(means) == 'new'] <- paste0("mean.", gsub("/", "_", group))
}
pdf( "histogram_all.pdf", width = 20, height = 9 )
#means = data.frame(target_id = abundances[["base/cdna"]]$target_id, mean.x10 = abundances[["base/cdna"]]$mean, mean.x20 = abundances[["base/illumina_probe"]]$mean);
merged <- merged.stack(means, id.vars = c("target_id"), var.stubs=c("mean"), sep=".")
setnames(merged, ".time_1", "group")
merged
wd <- "/home/umcg-jvanthuijl/umcg-jvanthuijl/source/R/gen_graphs/"
setwd(wd)
p1 <- ggplot(data=merged, aes(target_id, log10(mean))) +
geom_point(aes(color = (log10(mean) > threshold[2])))+
#geom_text(data = merged[log10(merged$mean) > threshold[2], ],
  #aes(label = target_id, x = target_id, y = log10(mean)), size = 4, hjust = -.1)+
scale_color_manual(values = c('black', 'red')) +
geom_hline(yintercept = threshold[2], color = "red", linetype = "dashed") + # threshold line
theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + #remove x axis
facet_grid( group ~. )
plot(p1)
#ggplot(merged, aes(x=group, y=log10(mean), fill=group)) + geom_boxplot()
#ggsave('plot_transcript.pdf', width = 10, height = 5, unit = 'in')
