## collateData.R by Rohan Maddamsetti.
## This script makes the final csv table that is used by analysis.R.

## import data.

thetaS.data <- read.csv("../data/thetaS.csv",header=T)
possible.synon.data <- read.csv("../data/possible_synonymous.csv", header=T)
ltee.synon.data <- read.csv("../data/ltee_synonymous_mutations.csv", header=T)
HGTs <- read.csv("../data/HGT_REL606.csv",header=T)

## change factors to character vectors (strings) so that they get merged
## properly.
thetaS.data$locus_tag <- sapply(thetaS.data$locus_tag,as.character)
possible.synon.data$locus_tag <- sapply(possible.synon.data$locus_tag,as.character)
possible.synon.data$gene <- sapply(possible.synon.data$gene,as.character)
ltee.synon.data$locus_tag <- sapply(ltee.synon.data$locus_tag,as.character)
ltee.synon.data$gene <- sapply(ltee.synon.data$gene,as.character)
HGTs$locus_tag <- sapply(HGTs$locus_tag,as.character)

## merge these data stepwise into a giant data frame.
murders <- merge(HGTs, thetaS.data, all=TRUE)
murders <- merge(murders, possible.synon.data, all=TRUE)
full.data <- merge(murders,ltee.synon.data,all=TRUE)

## write these data to file. PREPARE FOR ANALYSIS!
write.csv(full.data,file="../data/data.csv", row.names=F)
