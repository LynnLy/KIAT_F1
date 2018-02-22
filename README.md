# Reciprocal F1s in the KIAT project  

Supplement to Ruijuan's work:  

### Differential expression analysis using sleuth  
**F1_kallisto.sh** was used to trim and then pseudo-align the raw reads.  
**sleuth_F1.rmd** contains the sleuth object creation and analysis  

### Allele-specific expression analysis using MBASED  
**MBASED_F1.rmd** contains all preprocessing, filtering, and statistics   
**Run_MBASED.R** and **Run_MBASED_Beta.R** are the scripts for actually generating results from MBASED.  
The former assumes a binomial distribution of read counts and the latter assumes a beta-binomial distribution, with rho previously estimated in **MBASED_F1.rmd**  
**MBASED_F1_results.rmd** contains the results and analysis  
