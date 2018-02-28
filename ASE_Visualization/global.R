library(sleuth)
load("data/phasedData.Rdata")
load("data/aseGenesBeta.Rdata")
load("data/youngDEGenes_subset.Rdata")
so <- sleuth_load("data/SleuthObject")


visualData <- phasedData[, c("F1_414_RO", "F1_414_AO", "F1_415_RO", "F1_415_AO", "GeneID")]

visualDataTidy <- tidyr::gather(visualData, F1_414_RO, F1_414_AO, F1_415_RO, F1_415_AO, key = test, value = reads)
visualDataTidy <- tidyr::separate(visualDataTidy, col = test, into = c("discard", "cultivar", "allele"), 
                                  sep = "_", remove = TRUE, convert = FALSE)[, -2]

genes <- unique(visualDataTidy$GeneID)

DEGenes <- unique(youngDEGenes_subset$genes)
ASEDEGenes <- intersect(DEGenes, aseGenesBeta)
