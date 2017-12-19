#!/bin/env Rscript

#'
#' Option parsing
#'

usage <- function(mess = "", status = NULL) {
    message(mess)
    message("")
    message("Description:")
    message("This script is used to create a reference set from a directory with blck files.")
    message("")
    message("Usage:")
    message("Rscript nimbus_exomedepth_create_reference_set_from_directory.R --reference [reference_set file] --input [directory with block files]")
    message("")
    message("Options")
    message("--reference [file] the reference file to which to add the block file.")
    message("--input     [directory] the directory with block files, default: reference_set.txt")
    message("--help      []     prints this message")
    
    if(!is.null(status)){
        q(save="no", status=as.integer(status))
    }
}

# the option parsing library
library(getopt)

# Option parsing
spec <- matrix( c(
    'reference', 'r', 1, 'character',
    'input', 'i', 1, 'character',
    'help', 'h', 0, 'logical'
), byrow=TRUE, ncol=4 )
opt <- getopt(spec)

# Option checking
if( !is.null(opt$help)) { usage(status = 0) }

indir <- ifelse(!is.null(opt$input),{ opt$input }, {usage("No directory provided", status=2)})
fn_ref <- ifelse(!is.null(opt$reference),{ opt$reference }, {"reference_set.txt"})

ifelse(dir.exists(indir), NULL, {usage(sprintf("Input directory %s does not exist", indir), status=2)})

#'
#' Load essential libraries
#'

library(tidyverse)
library(stringr)

#'
#' Load the input files
#'

fn <- file.path(indir, list.files(indir, pattern = ".blck$"))
block <- bind_rows(lapply(fn, function(file) {
    read_tsv(file, col_names=c("key", "depth")) %>%
        mutate(sample = str_replace(str_replace(file, "^.*/", ""), "\\..*$", ""))
}))

#'
#' Generate the reference set
#'

reference <- block %>%
    spread(sample, depth)

#'
write_tsv(reference, fn_ref)