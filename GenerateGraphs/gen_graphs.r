library(ggplot2)
library(data.table)
library(splitstackshape)

rm(list=ls(all=TRUE))
wd <- "/home/umcg-jvanthuijl/umcg-jvanthuijl/kallisto"
setwd(wd)

# function to retrieve the abundance levels for all experiments within a directory dir.
calcAbundance <- function(dir, nrows = 20){
    outputAbundance = data.frame()
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

    # create a new dataframe, with the target_id and the mean of all samples
    outputAbundance <- data.frame(target_id = outputAbundance[,1],
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
corrr <- data.frame()

kmers=c("k11", "k21", "k31")
coverage=c("10x", "15x", "20x")
types=c("cdna",  "illumina_probe") #GENE REMOVED

samplegroups = paste(rep(types, each=6), rep(coverage, each=3), kmers, sep="/")
sampleprefixes = paste(rep(types, each=6), rep(coverage, each=3), kmers, sep="_")

#for (group in samplegroups){
#    abundances[[group]] = calcAbundance(paste(wd, "base", group, sep = "/"), 40000)
#    #   corrr$new <- abundances[[group]]$mean
#    # : names(corrr)[names(corrr) == "new"] <- group
#}
# read changed transcripts, in order to color their datapoints red, so it is
# known form the graph which points are expected to be outliers
changed_transcripts <- read.table("/home/umcg-jvanthuijl/umcg-jvanthuijl/results/simulated_reads/changed_transcripts.txt")

simulation_groups <- read.table("/home/umcg-jvanthuijl/umcg-jvanthuijl/results/simulated_reads/_diff_10/sim_rep_info.txt", header=T)

wd <- "/home/umcg-jvanthuijl/umcg-jvanthuijl/source/R/gen_graphs/"
setwd(wd)
diff_samplegroups <- paste("diff", samplegroups, sep="/")
colfunc <- colorRampPalette(c("red", "green"))
pdf( "diff_continuous.pdf", width = 10, height = 10 )

for (group in diff_samplegroups){
    diff1 <- calcDiffAbundance(paste(group, simulation_groups[simulation_groups$group==1,]$rep_id, sep = "/"), Inf)
    diff2 <- calcDiffAbundance(paste(group, simulation_groups[simulation_groups$group==2,]$rep_id, sep = "/"), Inf)
    plot(diff1$mean, diff2$mean, main = paste0(group, " n=10"), xlab = "base mean tpm", ylab = "diff mean tpm")
    # highlight the points  
    points(diff1$mean[as.numeric(rownames(changed_transcripts))],
           diff2$mean[as.numeric(rownames(changed_transcripts))],
           col=colfunc(13)[rank(changed_transcripts$V1)])

    subsetdiff1 <- data.frame(diff1[as.numeric(rownames(changed_transcripts)), ], expected = changed_transcripts$V1)
    subsetdiff1$sample = "Base"
    subsetdiff2 <- data.frame(diff2[as.numeric(rownames(changed_transcripts)), ], expected = changed_transcripts$V1)
    subsetdiff2$sample = "Diff"
    q <- rbind.data.frame(subsetdiff1, subsetdiff2)

    # ggplot(data = q, aes(target_id, log(mean), fill = sample, col = expected)) +
    #         geom_bar(stat = "identity", position="dodge") +
    #         coord_flip() + theme_bw() + labs( y = "Log mean tpm", x = "Transcript", title = group)

    # filename <- gsub("/", "-", group)
    # ggsave(paste0(filename, ".pdf"), device = "pdf", width = 10, height = 14)
}

stop()
M <- cor(corrr)
corrplot(M)

qplot(sd, data=abundances[["base/cdna/20x/k31"]], geom = "histogram", binwidth = 1)
plot(abundances[["base/cdna/10x/k31"]]$mean, abundances[["base/cdna/20x/k31"]]$mean)


# set lower and upper threshold for the graph.
# threshold = c(0.4, 1.3)
# 
# means = data.frame(target_id = abundances[["base/cdna"]]$target_id, mean.x10 = abundances[["base/cdna"]]$mean, mean.x20 = abundances[["base/illumina_probe"]]$mean);
# merged <- merged.stack(means, id.vars = c("target_id"), var.stubs=c("mean"), sep=".")
# setnames(merged, ".time_1", "group")
# merged
# 
# wd <- "/home/umcg-jvanthuijl/umcg-jvanthuijl/source/R/gen_graphs/"
# setwd(wd)
# p1 <- ggplot(data=merged, aes(target_id, log10(mean))) +
#   geom_point(aes(color = (log10(mean) > threshold[2])))+
#   geom_text(data = merged[log10(merged$mean) > threshold[2], ],
#             aes(label = target_id, x = target_id, y = log10(mean)), size = 4, hjust = -.1)+
#   scale_color_manual(values = c('black', 'red')) +
#   geom_hline(yintercept = threshold[2], color = "red", linetype = "dashed") + # threshold line
#   theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + #remove x axis
#   facet_grid( group ~. )
# plot(p1)
# 
# ggplot(merged, aes(x=group, y=log10(mean), fill=group)) + geom_boxplot()
# 
# ggsave('plot_transcript.pdf', width = 10, height = 5, unit = 'in')
