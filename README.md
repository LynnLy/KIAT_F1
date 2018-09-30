## Reciprocal F1 Observations

Script explanations below: 

## F1
### MBASED : Allele Specific Expression 
**MBASED_F1.rmd** contains all preprocessing, filtering, and statistics. Use this first to prepare data, and to calculate global reference bias, and to estimate `rho` for running MBASED using a beta-binomial distribution.    
**Run_MBASED.R** and **Run_MBASED_Beta.R** are the scripts for actually generating results from MBASED.  
The former assumes a binomial distribution of read counts and the latter assumes a beta-binomial distribution.  
**MBASED_F1_results.rmd** contains summary functions, results and analysis  

### kallisto : Read pseudoalignment
**F1_kallisto.sh** was used to trim and then pseudo-align the raw reads using `trimmomatic` and `kallisto`  
**subset_kallisto.R** contains a modified helper function `subset_kallisto()` from a defunct pull request for sleuth. This was needed to subset kallisto results, so that we could check for DE using ONLY genes that contain SNVs, for comparison to ASE results.  

### sleuth and edgeR : differential_expression  
`sleuth` uses `kallisto` bootstraps to check for differential expression while taking inferential variance (from the EM algorithm) into account. Our standard lab DE pipeline uses `edgeR` so I also included that.  
**sleuth_F1.rmd** contains `sleuth` creation and quality control. Use `sleuth_live(so)` to look at the results.    
**edgeR_F1.rmd** was used to check for DE using a second method, again with ONLY genes that contain SNVs  

### Shiny : ASE_Visualization
**ASE_Visualization** is a shiny app for visualizing the estimated kallisto read counts in each sample and the pooled allele specific reads for any given gene. You can choose which list of genes to look at: genes with ASE, DE, both, or all genes. 

## F2
**vcf_prep.rmd** Filtering and prepping VCF data for downstream QTL analysis  

### aseQTL : QTL using allele specific expression as the phenotype
**aseQTL_subset.rmd** uses 100 randomly sampled SNVs to do 1000 permutations on, for calculating LOD thresholds.  

To do: 
examine technical variance in DE genes called by edgeR but not sleuth (done but not reported yet)
