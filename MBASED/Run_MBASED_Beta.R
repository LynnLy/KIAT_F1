# R script for running the final analysis using a beta-binomial distribution.
# For analysis with the binomial distribution, use Run_MBASED.R
# Inputs: Pre-filtered, annotated, and phased data from MBASED_F1.rmd WITH dispersion parameter Rho
# The input data should be saved in "phasedData.Rdata"

setwd("~/repos/KIAT_F1/MBASED/Results_Beta")

rm(list = ls())

library(MBASED)
library(tidyverse)

TwoSample <- function(annotatedData, mySNVs, genotype1, genotype2, numSim = 0){
  RO1 <- paste(genotype1, "RO", sep = "_")
  AO1 <- paste(genotype1, "AO", sep = "_")
  RO2 <- paste(genotype2, "RO", sep = "_")
  AO2 <- paste(genotype2, "AO", sep = "_")
  RAB1 <- paste(genotype1, "refBias", sep = "_")
  RAB2 <- paste(genotype2, "refBias", sep = "_")
  DISP1 <- paste(genotype1, "disp", sep = "_")
  DISP2 <- paste(genotype2, "disp", sep = "_")
  
  mySample <- SummarizedExperiment(
    assays = list(lociAllele1Counts = matrix(c(annotatedData[, RO1], annotatedData[, RO2]), ncol = 2,
                                             dimnames = list(names(mySNVs), c(genotype1, genotype2))),
                  
                  lociAllele2Counts = matrix(c(annotatedData[, AO1], annotatedData[, AO2]), ncol = 2,
                                             dimnames = list(names(mySNVs), c(genotype1, genotype2))),
                  
                  lociAllele1CountsNoASEProbs = matrix(c(annotatedData[, RAB1], annotatedData[, RAB2]),
                                                       ncol=2, dimnames=list(names(mySNVs), c(genotype1, genotype2))),
                  lociCountsDispersions = matrix(c(annotatedData[, DISP1], annotatedData[, DISP2]), 
                                                       ncol=2, dimnames=list(names(mySNVs), c(genotype1, genotype2)))
                  ),
                  
    rowRanges=mySNVs)
  
  MBASEDOutput <- runMBASED(
    ASESummarizedExperiment = mySample,
    isPhased = TRUE,
    numSim = numSim,
    BPPARAM = SnowParam(workers = 8) # Default: No paralellization
  )
  
  return(MBASEDOutput)
} 

runOnSubset <- function(annotatedData, index){
  annotatedData.trimmed <- annotatedData[index, ]
  
  mySNVs.trimmed <- GRanges(
    seqnames = annotatedData.trimmed$CHROM,
    ranges = IRanges(start = annotatedData.trimmed$POS, width = 1),
    aseID = as.vector(annotatedData.trimmed$GeneID),
    allele1 = annotatedData.trimmed$REF,
    allele2 = annotatedData.trimmed$ALT)
  
  return(TwoSample(annotatedData.trimmed, mySNVs.trimmed, "F1_414", "F1_415", numSim = 1000000))
}

# IMPORTANT: 2-SAMPLE ANALYSIS IS NOT SYMMETRIC. This function is for using F1_415 as sample 1. 
runOnSubsetReverse <- function(annotatedData, index){
  annotatedData.trimmed <- annotatedData[index, ]
  
  mySNVs.trimmed <- GRanges(
    seqnames = annotatedData.trimmed$CHROM,
    ranges = IRanges(start = annotatedData.trimmed$POS, width = 1),
    aseID = as.vector(annotatedData.trimmed$GeneID),
    allele1 = annotatedData.trimmed$REF, 
    allele2 = annotatedData.trimmed$ALT)
  
  return(TwoSample(annotatedData.trimmed, mySNVs.trimmed, "F1_415", "F1_414", numSim = 1000000))
}

load("../phasedData.Rdata")
phasedData <- phasedData[, -c(5:18)]

SNVsPerGene <- phasedData %>% group_by(GeneID) %>% tally()

# Gather all the genes with over 20 SNVs, as we don't have enough memory on Whitney to compute them.
over20SNVs <- SNVsPerGene[SNVsPerGene$n > 20, ] # There are 33 of them. We can return to them later. 

phasedData <- phasedData[!(phasedData$GeneID %in% over20SNVs$GeneID), ]


FindIndexes <- function(data, start, size, end = nrow(data)) {
  # Required so that we don't separate SNVs that are from the same gene when subsetting our data
  # Start = start index
  # Size = Desired approximate interval length
  i <- start + size
  
  if(i >= end) {return(end)} # Reached the end of the dataset
  
  while(data$newGene[i] == FALSE) {i <- i + 1}
  
  return(i - 1)
}

startIndexes <- rep(1,42)
endIndexes <- rep(1, 42)
i = 1
while(i < 42) {
  endIndexes[i] <- FindIndexes(phasedData, start = startIndexes[i], size = 1000)
  i <- i + 1
  startIndexes[i] <- endIndexes[i-1] + 1
}

for (i in 1:41) {
  old <- Sys.time()
  start <- startIndexes[i]
  end <- endIndexes[i]
  
  #subsetName <- paste0("BETA.MBASED.F1.414.vs.F1.415.", start, "to", end)
  subsetName <- paste0("BETA.MBASED.F1.415.vs.F1.414.", start, "to", end)
  print(paste0(i, "th Index. Working on ", subsetName))
  #assign(subsetName, runOnSubset(phasedData, index = start:end))
  assign(subsetName, runOnSubsetReverse(phasedData, index = start:end))
  save(list = c(subsetName), file = paste0(subsetName, ".Rdata"))
  rm(list = c(subsetName))
  
  new <- Sys.time() - old 
  print(paste0("Time elapsed: ", new))
}
