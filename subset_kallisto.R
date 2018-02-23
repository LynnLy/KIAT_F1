# Helper function to subset a kallisto

# Taken from https://github.com/pachterlab/sleuth/pull/150/files

filter_bootstraps <- function(bs, target_ids) { 
  # I made this function because the bootstrap filter in below was not working
  filteredBootstraps <- bs[[1]][which(bs[[1]]$target_id %in% target_ids), ]
  for(i in 2:length(bs)) {
    filteredBootstraps <- rbind(filteredBootstraps, bs[[i]][which(bs[[i]]$target_id %in% target_ids), ])
  }
  return(filteredBootstraps)
}

testboot <- filter_bootstraps(testKal$bootstrap, geneIDs)

subset_kallisto_custom <- function(obj, target_ids) {
  stopifnot(is(target_ids, "character"))
  stopifnot(all(target_ids %in% obj$abundance$target_id))

  subset_num <- length(unique(target_ids))
  new_obj <- obj
  new_obj$abundance <- new_obj$abundance[which(new_obj$abundance$target_id %in% target_ids), ] # I changed this line
  #new_obj$bootstrap <- lapply(new_obj$bootstrap, function(bs) {
  # bs[which(bs$target_id %in% target_ids), ]
  #})
  new_obj$bootstrap <- filter_bootstraps(new_obj$bootstrap, target_ids)
  
  
  attr(new_obj, "num_targets") <- subset_num
  excluded_ids <- obj$abundance$target_id[which(!(obj$abundance$target_id %in% target_ids))]
  if(length(new_obj$excluded_ids) == 0) {
    new_obj$excluded_ids <- excluded_ids
  } else {
    new_obj$excluded_ids <- c(new_obj$excluded_ids, excluded_ids)
  }
  new_obj
}

#' Write a kallisto object to HDF5
#'
#' Write a kallisto object to HDF5.
#'
#' @param kal the kallisto object to write out
#' @param fname the file name to write out to
#' @param overwrite whether the file should be overwritten if it exists
#' @param write_bootstrap whether the bootstraps should be written to file
#' @param compression an integer between 0 and 7 that indicates the level of compression
#'    to use, with 0 being no compression and 7 being the highest supported by this method.
#'    The default of 6 is a good choice for most applications.
#' @return the kallisto object \code{kal} invisibly.
#' @importFrom rhdf5 h5write.default
#' @importFrom rhdf5 h5write
#' @export
write_kallisto_hdf5 <- function(kal, fname, overwrite = TRUE,
                                write_bootstrap = TRUE, compression = 6L) {
  stopifnot(is(kal, "kallisto"))
  stopifnot(is(fname, "character"))
  stopifnot(is(compression, "integer"))
  stopifnot(length(compression) == 1)
  
  # TODO: ensure that all bootstraps are sorted according to abundance
  if (compression < 0 || compression > 7 ) {
    stop("'compression' must be an integer between 0 and 7")
  }
  
  fname <- path.expand(fname)
  
  if (file.exists(fname)) {
    if (overwrite) {
      warning(paste0("'", fname, "' already exists. Overwriting."))
      file.remove(fname)
    } else {
      stop(paste0("'", fname, "' already exists."))
    }
  }
  
  if (!rhdf5::h5createFile(fname)) {
    stop(paste0("Error: Couldn't open '", fname, "' to write out."))
  }
  
  dims <- c(nrow(kal$abundance), 1)
  cat("dims: ", dims, "\n")
  
  # write out auxilary info
  rhdf5::h5createGroup(fname, "aux")
  
  stopifnot(rhdf5::h5createDataset(fname, "aux/ids", dims = dims,
                                   storage.mode = "character", size = 100,
                                   level = compression))
  
  if (write_bootstrap) {
    rhdf5::h5write(length(kal$bootstrap), fname, "aux/num_bootstrap")
  } else {
    rhdf5::h5write(0L, fname, "aux/num_bootstrap")
  }
  
  rhdf5::h5write(kal$abundance$target_id, fname, "aux/ids")
  rhdf5::h5write(kal$abundance$eff_len, fname, "aux/eff_lengths")
  rhdf5::h5write(kal$abundance$len, fname, "aux/lengths")
  rhdf5::h5write(kal$fld, fname, "aux/fld")
  rhdf5::h5write(kal$bias_normalized, fname, "aux/bias_normalized")
  rhdf5::h5write(kal$bias_observed, fname, "aux/bias_observed")
  rhdf5::h5write(attributes(kal)$num_processed, fname, "aux/num_processed")
  rhdf5::h5write(attributes(kal)$index_version, fname, "aux/index_version")
  rhdf5::h5write(attributes(kal)$kallisto_version, fname, "aux/kallisto_version")
  rhdf5::h5write(attributes(kal)$start_time, fname, "aux/start_time")
  rhdf5::h5write(kal$abundance$est_counts, fname, "est_counts")
  
  if (write_bootstrap && length(kal$bootstrap) > 0) {
    rhdf5::h5createGroup(fname, "bootstrap")
    for (i in seq_along(kal$bootstrap)) {
      bs <- kal$bootstrap[[i]]$est_counts
      bs_name <- paste0("bootstrap/bs", i-1)
      rhdf5::h5write(bs, fname, bs_name)
    }
  }
  
  invisible(kal)
}