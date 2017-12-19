#!/bin/env Rscript


usage <- function(mess = "", status = NULL) {
    message(mess)
    message("")
    message("Description:")
    message("This script is used to add a single block file to a reference set.")
    message("")
    message("Usage:")
    message("Rscript nimbus_exomedepth.R --reference [reference_set file] --input [block file]")
    message("")
    message("Options")
    message("--reference [file] the reference file to which to add the block file")
    message("--input     [file] the block file to add the reference")
    message("--help      []     prints this message")
    
    if(!is.null(status)){
        q(save="no", status=as.integer(status))
    }
}

#'
#' Option parsing
#'

# the option parsing library
library(getopt)

spec <- matrix( c(
    'reference', 'r', 1, 'character',
    'input', 'i', 1, 'character',
    'threshold', 't', 1, 'integer',
    'output', 'o', 1, 'character',
    'bed', 'b', 1, 'character',
    'help', 'h', 0, 'logical'
), byrow=TRUE, ncol=4 )
opt <- getopt(spec)

# Option checking
if( !is.null(opt$help)) { usage(status = 0) }
fn_tst <- ifelse(!is.null(opt$input),{ opt$input }, {usage("No input file provided", status=2)})
fn_ref <- ifelse(!is.null(opt$reference),{ opt$reference }, {usage("No reference file specified", status=2)})
ifelse(file.exists(fn_tst), NULL, {usage(sprintf("Input file %s does not exist", fn_tst), status=2)})

#'
#' Load the required libraries
#'

library(tidyverse)
library(stringr)

#'
#' Load the data
#'

#'
block <- read_tsv(fn_tst, col_names=c("key", str_replace(str_replace(fn_tst, "^.*/", ""), "\\..*$", "")))

#'
#' Add the block file to the reference if it exists
#' 

if(!file.exists(fn_ref)){
    reference <- read_tsv(fn_ref)
    block <- left_join(reference, block, by = c("key"))
}

#' write the reference
write_tsv(block, fn_ref)