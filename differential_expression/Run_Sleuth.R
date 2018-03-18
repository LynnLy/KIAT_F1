# This script is to run sleuth using multiple cores, since you can only use one core when not using the command line

library(sleuth)

# set the number of available cores to 4
options(mc.cores = 4L)

# collect the file paths
sample_id <- dir(file.path("kallisto/results"))
kal_dirs <- file.path("kallisto/results", sample_id, "kallisto") # This line was used for sleuth with all the genes
# kal_dirs <- file.path("kallisto/results_subset", sample_id, "kallisto") # This line was used for sleuth with a subset of the genes
kal_dirs

# collect the metadata and combine it with file paths
metadata <- read.csv(file.path("F1_summary.csv"), header = TRUE, stringsAsFactors=FALSE)
metadata <- metadata[1:6, ]# young tissue only
metadata <- dplyr::mutate(metadata, path = kal_dirs)
metadata <- dplyr::select(metadata[1:6, ], sample = Sample.ID, cultivar, path)

# filtering step
custom_filter <- function(row, min_reads = 5, min_prop = 0.47) {
  (mean(row >= min_reads) >= min_prop) 
}

so <- sleuth_prep(metadata, extra_bootstrap_summary = TRUE, filter = custom_filter)

so <- sleuth_fit(so, ~cultivar, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

so <- sleuth_wt(so, 'cultivar415F1')

sleuth_save(so, file = "sleuthResults.Rdata")