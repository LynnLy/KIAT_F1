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
    numSim = numSim
    #BPPARAM = MulticoreParam(workers = 9) # Default: No paralellization
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

load("phasedData.Rdata")

BMBASED.F1.414.vs.F1.415A <- runOnSubset(phasedData, 1:4997)
save(BMBASED.F1.414.vs.F1.415A, file = "BMBASED.F1.414.vs.F1.415A.Rdata")

BMBASED.F1.414.vs.F1.415B1 <- runOnSubset(phasedData, 4998:7501)
save(BMBASED.F1.414.vs.F1.415B1, file = "BMBASED.F1.414.vs.F1.415B1.Rdata")

BMBASED.F1.414.vs.F1.415B2 <- runOnSubset(phasedData, 7502:10002)
save(BMBASED.F1.414.vs.F1.415B2, file = "BMBASED.F1.414.vs.F1.415B2.Rdata")

BMBASED.F1.414.vs.F1.415C1 <- runOnSubset(phasedData, 10003:12501)
save(BMBASED.F1.414.vs.F1.415C1, file = "BMBASED.F1.414.vs.F1.415C1.Rdata")

BMBASED.F1.414.vs.F1.415C2 <- runOnSubset(phasedData, 12502:15006)
save(BMBASED.F1.414.vs.F1.415C2, file = "BMBASED.F1.414.vs.F1.415C2.Rdata")

print("Working on E1")

BMBASED.F1.414.vs.F1.415E1 <- runOnSubset(phasedData, 20005:22507)
save(BMBASED.F1.414.vs.F1.415E1, file = "BMBASED.F1.414.vs.F1.415E1.Rdata")
print("working on F2")

BMBASED.F1.414.vs.F1.415F2 <- runOnSubset(phasedData, 27505:30004)
save(BMBASED.F1.414.vs.F1.415F2, file = "BMBASED.F1.414.vs.F1.415F2.Rdata")

print("working on H")

BMBASED.F1.414.vs.F1.415H <- runOnSubset(phasedData, 35004:37502)
save(BMBASED.F1.414.vs.F1.415H, file = "BMBASED.F1.414.vs.F1.415H.Rdata")

print("working on I")

BMBASED.F1.414.vs.F1.415I <- runOnSubset(phasedData, 37503:41404)
save(BMBASED.F1.414.vs.F1.415I, file = "BMBASED.F1.414.vs.F1.415I.Rdata")

print("working on E2")

BMBASED.F1.414.vs.F1.415E2 <- runOnSubset(phasedData, 22508:24002)
save(BMBASED.F1.414.vs.F1.415E2, file = "BMBASED.F1.414.vs.F1.415E2.Rdata")

print("working on E3")

BMBASED.F1.414.vs.F1.415E3 <- runOnSubset(phasedData, 24003:25001)
save(BMBASED.F1.414.vs.F1.415E3, file = "BMBASED.F1.414.vs.F1.415E3.Rdata")

print("working on F1") 

BMBASED.F1.414.vs.F1.415F1 <- runOnSubset(phasedData, 25002:27504)
save(BMBASED.F1.414.vs.F1.415F1, file = "BMBASED.F1.414.vs.F1.415F1.Rdata")

