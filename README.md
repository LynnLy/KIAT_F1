# Reciprocal F1s in the KIAT project  

To do: generalize the kallisto-sleuth workflow using Snakemake  
QTL using ASE as the trait of interest ('aseQTL')  
examine technical variance in DE genes called by edgeR but not sleuth  

Supplement to Ruijuan's work:  

### Allele-specific expression analysis using MBASED  
**MBASED_F1.rmd** contains all preprocessing, filtering, and statistics. Use this first to prepare data, and to calculate global reference bias, and to estimate `rho` for running MBASED using a beta-binomial distribution.    
**Run_MBASED.R** and **Run_MBASED_Beta.R** are the scripts for actually generating results from MBASED.  
The former assumes a binomial distribution of read counts and the latter assumes a beta-binomial distribution.  
**MBASED_F1_results.rmd** contains the results and analysis  

### Visualization using Shiny
**ASE_Visualization** is a shiny app for visualizing the estimated kallisto read counts in each sample and the pooled allele specific reads for any given gene. You can choose which list of genes to look at: genes with ASE, DE, both, or all genes.   

### Differential expression analysis using sleuth  
**F1_kallisto.sh** was used to trim and then pseudo-align the raw reads using `trimmomatic` and `kallisto`  
**sleuth_F1.rmd** contains `sleuth` creation and QC  
**subset_kallisto.R** contains a modified helper function `subset_kallisto()` from a defunct pull request for sleuth. This was needed to subset kallisto results, so that we could check for DE using ONLY genes that contain SNVs, for comparison to ASE results.  

#### Using edgeR  
**edgeR_F1.rmd** was used to check for DE using a second method, again with ONLY genes that contain SNVs  
