## analysis.R by Rohan Maddamsetti.

## This R script analyzes data.csv, a file containing data on the number
## of HGT sequences from Escherichia or Shigella in the hgt.fst data set
##(Smillie 2011) that BLAST onto E. coli B str. REL606 with an E-value smaller
## than 10E-10, as well as thetaS values calculated using omegamap
## (Wilson 2006; Martincorena 2012) on a set of
## panortholog (Cooper 2010) alignments from 10 E. coli genomes:
## E-coli-B-REL606, E-coli-O17-K52-H18-UMN026, E-coli-O26-H11-11368,
## E-coli-O111-H-11128, E-coli-O9-HS, E-coli-O157-H7-Sakai-EHEC,
## E-coli-SECEC-SMS-3-5, E-coli-BL21-DE3, E-coli-O6-K15-H31-536-UPEC,
## and E-coli-K-12-MG1655. data.csv also contains genomic data on
## synonymous substitutions in all long term lines
## after 40,000 generations had elapsed, as well as the number of all possible
##kinds of synonymous point mutations in every gene in E. coli B str. REL606.


addTotalMutationsColumn <- function(the.data) {
  ## Add a column to the data frame counting the total number of synonymous
  ## substitutions in the LTEE at 40K at each locus.
  ## The lambda function turns NA values to zeros so the summation works nicely.
  ## This function is called in prepareData(filename).
  data <- the.data
  data$total.mutations <- sapply(data$G.to.A,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$G.to.T,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$G.to.C,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$A.to.G,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$A.to.T,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$A.to.C,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$T.to.G,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$T.to.A,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$T.to.C,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$C.to.G,function (x) ifelse(is.na(x),0,x)) +
  sapply(data$C.to.A,function (x) ifelse(is.na(x),0,x)) +
    sapply(data$C.to.T,function (x) ifelse(is.na(x),0,x))

  return(data)
}

prepareData <- function(filename) {
  ## This function reads in the data from data.csv, and prepares it for
  ## analysis.

  ##import data.

  the.data <- read.csv(filename,header=T)

  ## add the total mutations column.
  the.data <- addTotalMutationsColumn(the.data)

  ##remove entries with missing thetaS values.
  the.data <- the.data[!is.na(the.data$thetaS),]

##remove entries with missing "possible.G.to.A" etc. values--
##these are either pseudogenes, or genes containing selenocysteine.

  the.data <- the.data[!is.na(the.data$possible.G.to.A),]

##sort the data frame by thetaS.

  the.data <- the.data[with(the.data,order(thetaS)),]

## HGT.hits are converted to TRUE/FALSE values, since the number of hits is not
## a reliable measure of the true number of HGT events including this sequence
## (personal communication with Chris Smillie).
  the.data$HGT.hits <- sapply(the.data$HGT.hits, function (x) ifelse(is.na(x),FALSE,TRUE))
  
  return(the.data)
  
}

plotMutationSpectra <- function (the.data) {
## This function plots the mutation spectra of the 12 LTEE clones.
  
  library(ggplot2)

  genomes <- levels(the.data$genome)
 
  spectrum.data <- data.frame(Genome=genomes, AT.to.TA=rep(NA,12),
                              AT.to.CG=rep(NA,12),
                              AT.to.GC=rep(NA,12),
                              CG.to.TA=rep(NA,12),
                              CG.to.AT=rep(NA,12),
                              CG.to.GC=rep(NA,12),
                              total.mutations=rep(NA,12))

  for (genome in genomes) {
    total.AT.to.TA <- sum(the.data[the.data$genome==genome,]$A.to.T,na.rm=TRUE) + sum(the.data[the.data$genome==genome,]$T.to.A,na.rm=TRUE)
    
    total.AT.to.CG <- sum(the.data[the.data$genome==genome,]$A.to.C,na.rm=TRUE) + sum(the.data[the.data$genome==genome,]$T.to.G,na.rm=TRUE)
    
    total.AT.to.GC <- sum(the.data[the.data$genome==genome,]$A.to.G,na.rm=TRUE) + sum(the.data[the.data$genome==genome,]$T.to.C,na.rm=TRUE)
    
    total.CG.to.TA <- sum(the.data[the.data$genome==genome,]$C.to.T,na.rm=TRUE) + sum(the.data[the.data$genome==genome,]$G.to.A,na.rm=TRUE)
    
    total.CG.to.AT <- sum(the.data[the.data$genome==genome,]$C.to.A,na.rm=TRUE) + sum(the.data[the.data$genome==genome,]$G.to.T,na.rm=TRUE)
    
    total.CG.to.GC <- sum(the.data[the.data$genome==genome,]$C.to.G,na.rm=TRUE) + sum(the.data[the.data$genome==genome,]$G.to.C,na.rm=TRUE)

    spectrum.data[spectrum.data$Genome==genome,]$total.mutations <- sum(total.AT.to.TA, total.AT.to.CG, total.AT.to.GC, total.CG.to.TA, total.CG.to.AT, total.CG.to.GC)
    spectrum.data[spectrum.data$Genome==genome,]$AT.to.TA <- total.AT.to.TA
    spectrum.data[spectrum.data$Genome==genome,]$AT.to.CG <- total.AT.to.CG
    spectrum.data[spectrum.data$Genome==genome,]$AT.to.GC <- total.AT.to.GC
    spectrum.data[spectrum.data$Genome==genome,]$CG.to.TA <- total.CG.to.TA
    spectrum.data[spectrum.data$Genome==genome,]$CG.to.AT <- total.CG.to.AT
    spectrum.data[spectrum.data$Genome==genome,]$CG.to.GC <- total.CG.to.GC
  }

  total <- sum(spectrum.data$total.mutations)
  rows <- c("AT.to.TA", "AT.to.CG", "AT.to.GC", "CG.to.TA", "CG.to.AT", "CG.to.GC")

  spectrum.barplot.data <- data.frame(Genome=c(),Mutation=c())
  for (g in genomes) {
    for (j in rows) {
      reps <- subset(spectrum.data[spectrum.data$Genome==g,], select=c(j))
      if (reps > 0) {
        current <- data.frame(Genome=rep(g,reps), Mutation=rep(j, reps))
        spectrum.barplot.data <- rbind(spectrum.barplot.data,current)
      }
    }
  }
  spectrum.barplot.data$Genome <- sapply(spectrum.barplot.data$Genome, function (x) substr(x,1,5))
  spectrum.plot <- ggplot(spectrum.barplot.data, aes(x=Genome)) + geom_bar(aes(fill=Mutation))
  spectrum.plot
  ggsave("../results/mutation_spectrum.pdf")
}

ks.analysis <- function (the.data) {
  ## For each set of data (all data, non-mutators, MMR mutators, mutT mutators)
  ## do the following: 1) make a uniform cdf on mutation rate per base.
  ## 2) make a thetaS cdf. 3) make an empirical cdf of mutations per gene.
  ## do K-S tests for goodness of fit of the empirical cdf with the cdfs for
  ## the uniform cdf and thetaS cdf hypotheses.
  
  relevant.data <- subset(the.data, select=c("locus_tag","gene", "thetaS", "total.mutations", "gene.length", "HGT.hits"))
  unique.relevant.data <- unique(relevant.data)

  genome.length <- sum(unique.relevant.data$gene.length)
  
  ## Calculate the empirical distribution of synonymous substitutions per gene.
  mutation.total <- sum(unique.relevant.data$total.mutations)
  empirical.cdf <- cumsum(unique.relevant.data$total.mutations/mutation.total)
  ## Null hypothesis: probability of a mutation per base is uniform.
  null.cdf <- cumsum(unique.relevant.data$gene.length/genome.length)
  ## Alternative hypothesis: mutation rate is proportional to thetaS.
  ## alternative 1: thetaS is the mutation rate per base pair.
    thetaS.is.per.bp.total <- sum(unique.relevant.data$thetaS*unique.relevant.data$gene.length)
    alt1.cdf <- cumsum(unique.relevant.data$thetaS*unique.relevant.data$gene.length/thetaS.is.per.bp.total)
  ## alternative 2: thetaS is the mutation rate per gene.
  thetaS.is.per.gene.total <- sum(unique.relevant.data$thetaS)
  alt2.cdf <- cumsum(unique.relevant.data$thetaS/thetaS.is.per.gene.total)
  ## Do Kolmogorov-Smirnov tests for goodness of fit.
  print(ks.test(empirical.cdf, null.cdf, simulate.p.value=TRUE))
  print(ks.test(empirical.cdf,alt1.cdf, simulate.p.value=TRUE))
  print(ks.test(empirical.cdf,alt2.cdf, simulate.p.value=TRUE))

  results.to.plot <- data.frame(locus_tag=unique.relevant.data$locus_tag, gene=unique.relevant.data$gene, thetaS=unique.relevant.data$thetaS, empirical=empirical.cdf,null=null.cdf,alt1=alt1.cdf,alt2=alt2.cdf, HGT.hits=unique.relevant.data$HGT.hits)

  return(results.to.plot)
}

makeFigure <- function(the.results.to.plot) {
## This function generates the first panel that I want
## to use to illustrate how there is an association with thetaS with HGT,
## but not with synonymous substitution rates in the LTEE.
## Some post-processing with Illustrator will be needed (basically to add
## thetaS as a greek symbol on the axes).

  library(ggplot2)
  
  ## Plot thetaS across all loci.
  ## On this same plot, plot HGT.hits as red dots.
  ## '0' in the axes is a placeholder for adding the greek letter theta in
  ## post-processing with Illustrator.
  
  ## for plotting convienence, add an index to the data frame.
  the.results.to.plot$index <- 1:length(the.results.to.plot$gene)
  
  plot <- ggplot(the.results.to.plot, aes(x=index)) + geom_point(aes(y=empirical, group=HGT.hits, colour=HGT.hits, size=3)) + geom_line(aes(y=null, linetype=2)) + geom_line(aes(y=alt1, linetype=4,size=2)) + scale_colour_manual(name="HGT hit", values=c("light green", "black")) + scale_x_continuous('Genes ranked by 0',limits=c(0,2900)) + scale_y_continuous('Cumulative proportion of synonymous mutations',limits=c(0,1)) + opts(panel.grid.minor = theme_blank(), panel.background = theme_rect(colour = "white"), axis.line = theme_segment(), legend.position="none", axis.text.y=theme_text(colour="black", size=16), axis.text.x=theme_text(colour="black", size=16), axis.title.x=theme_text(size=20), axis.title.y=theme_text(angle=90,size=20))
  
  plot
  ggsave("../results/figure.pdf")
  ggsave("../results/figure.ps")
  
}

## Import the data, and get it ready for analysis.
data <- prepareData("../data/data.csv")

## I classify the different kinds of mutators based on the mutational spectrum:
## Ara+6 for example has mutations in mutL, mutS, and mutT, but acts as a mutT
## mutator. Running the function plotMutationSpectra shows this nicely.

plotMutationSpectra(data)

MMR.genomes <- c("Ara-3_REL10979", "Ara+3_REL10953", "Ara-4_REL10944", "Ara-2_REL11036")
mutT.genomes <- c("Ara-1_REL10938", "Ara+6_REL10985")
nonmutator.genomes <- c("Ara-5_REL10947", "Ara-6_REL11005", "Ara+1_REL11008", "Ara+2_REL10950", "Ara+4_REL10956", "Ara+5_REL10982")
  
MMR.data <- data[data$genome %in% MMR.genomes,]
mutT.data <- data[data$genome %in% mutT.genomes,]
nonmutator.data <- data[data$genome %in% nonmutator.genomes,]

mmr.results <- ks.analysis(MMR.data)
mutT.results <- ks.analysis(mutT.data)
nonmutator.results <- ks.analysis(nonmutator.data)
results <- ks.analysis(data)

makeFigure(results)

## Do analysis of HGT data.
## Do a Mann-Whitney U-test to see if the distribution of thetaS values among
## loci hit by HGT is different from the distribution of thetaS values among
## loci not hit by HGT.

relevant.data <- subset(data, select=c("locus_tag","HGT.hits","thetaS"))

## remove duplicate entries.
unique.relevant.data <- unique(relevant.data)

HGT.trues <- unique.relevant.data[unique.relevant.data$HGT.hits==TRUE,]$thetaS
HGT.falses <- unique.relevant.data[unique.relevant.data$HGT.hits==FALSE,]$thetaS
wilcox.test(HGT.trues,HGT.falses,paired=FALSE,alternative=c("greater"))

## repeat HGT analysis with Martincorena's thetaS estimates.
martincorena.data <- read.csv("../data/Martincorena_Maddamsetti_thetaS_estimates.csv", header=T)
martincorena.HGT <- merge(martincorena.data, unique.relevant.data, by=("locus_tag")
martincorena.HGT.trues <- martincorena.HGT[martincorena.HGT$HGT.hits==TRUE,]$Martincorena_thetaS)
martincorena.HGT.trues <- martincorena.HGT.trues[!is.na(martincorena.HGT.trues)]
martincorena.HGT.falses <- martincorena.HGT[martincorena.HGT$HGT.hits==FALSE,]$Martincorena_thetaS
martincorena.HGT.falses <- martincorena.HGT.falses[!is.na(martincorena.HGT.falses)]
wilcox.test(martincorena.HGT.trues,martincorena.HGT.falses,paired=FALSE,alternative=c("greater"))
