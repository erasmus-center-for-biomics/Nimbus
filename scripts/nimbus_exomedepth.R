#!/bin/env Rscript

#'
#' Option parsing
#'

usage <- function(mess = "", status = NULL) {
    message(mess)
    message("")
    message("Description:")
    message("This script runs ExomeDepth* for Nimbus blck files.")
    message("")
    message("*) Plagnol,V. et al. (2012) A robust model for read count data in exome sequencing experiments and implications for copy number variant calling. Bioinformatics, 28, 2747–54.")
    message("")
    message("Usage:")
    message("Rscript nimbus_exomedepth.R --reference [reference_set file] --input [block file] --threshold [minimum number of reads per amplicon] --output [exome depth output file]")
    message("")
    message("Options")
    message("--reference  [file] the reference file")
    message("--input      [file] the input block file")
    message("--output     [file] the output exome depth file, default: exomedepth_cnv_calls.txt")
    message("--threshold  [int]  the coverage threshold for amplicons in the reference set to be included, default: 10")
    message("--nbins      [int]  the number of bins for the selectt reference set, default: 10000")
    message("--transition [double]  the transition parameter for ExomeDepth, default: 10^-4")
    message("--selection  [file] The output for the dataset selection, default: selection_statistics.txt")
    message("--amplicons  [file] The amplicon fold changes and z-scores compared to selected references, default: amplicon_statistics.txt")
    message("--help       []     prints this message")
    
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
    'threshold', 't', 2, 'integer',
    'transition', 'x', 2, 'float',
    'nbins', 'y', 2, 'integer',
    'selection', 'z', 2, 'character',
    'amplicon', 'a', 2, 'character',
    'output', 'o', 2, 'character',
    'help', 'h', 0, 'logical'
), byrow=TRUE, ncol=4 )
opt <- getopt(spec)

# Option checking
if( !is.null(opt$help)) { usage(status = 0) }

fn_tst <- ifelse(!is.null(opt$input),{ opt$input }, {usage("No input file provided", status=2)})
fn_ref <- ifelse(!is.null(opt$reference),{ opt$reference }, {usage("No reference file specified", status=2)})
fn_out <- ifelse(!is.null(opt$output),{ opt$output }, {"exomedepth_cnv_calls.txt"})
transition <- ifelse(!is.null(opt$transition),{ opt$transition }, {10^-4})
nbins <- ifelse(!is.null(opt$nbins),{ opt$nbins }, {10000})
fn_selection <- ifelse(!is.null(opt$selection),{ opt$selection }, {"selection_statistics.txt"})
fn_amplicon <- ifelse(!is.null(opt$amplicon),{ opt$amplicon }, {"amplicon_statistics.txt"})
read_thr <- ifelse(!is.null(opt$threshold),{ opt$threshold }, 10)

#' Check the input files
ifelse(file.exists(fn_tst), NULL, {usage(sprintf("Input file %s does not exist", fn_tst), status=2)})
ifelse(file.exists(fn_ref), NULL, {usage(sprintf("Reference set %s does not exist", fn_ref), status=2)})

#'
#' Load the required libraries
#'

library(IRanges)
library(ExomeDepth)
library(tidyverse)
library(stringr)

#'
#' Load the test file
#'

block <- read_tsv(fn_tst, col_names=c("key", "depth"))

#'
#' Load the reference file
#'

reference <- read_tsv(fn_ref)

#' filter the reference
tf <- apply(reference[,2:ncol(reference)], 1, function(x) all(x > read_thr))
reference <- reference[tf,]

#' filter the new data
idx <- match(reference$key, block$key)
block <- block[idx,]

#'
#' Select the most appropriate sets
#'

choice <-  ExomeDepth::select.reference.set(
    test.counts      = block$depth,
    reference.counts = as.matrix(reference[,2:ncol(reference)]),
    n.bins.reduced   = nbins)

selection <- as_tibble(choice$summary.stats) %>%
    mutate(`in reference set` = ref.samples %in% choice$reference.choice)

write_tsv(selection, fn_selection)

#'
#' Run exome depth
#'

v_reference <- rowSums(as.matrix(reference[,v$reference.choice]))
object <- new(
    'ExomeDepth', 
    test = block$depth, 
    reference = v_reference, 
    formula = 'cbind(test, reference) ~ 1' )

# get the chromosome/start/end/strand
matcher <- str_match(block$key, "^(.+):([0-9]+)-([0-9]+)\\(([+,-])\\)")

#'
#' Call the CNVs from the data
#'

result <- ExomeDepth::CallCNVs(
    x = object,
    transition.probability = transition,
    chromosome = matcher[,2],
    start = parse_integer(matcher[,3]) + 1,
    end = parse_integer(matcher[,3]) + 2,
    name = block$key)

cnv_calls <- as_tibble(result@CNV.calls) %>%
    arrange(desc(BF)) %>%
    mutate(id = str_replace(id, "^chr", ""))

# write the output
write_tsv(cnv_calls, fn_out)

#'
#' Determine the fold-changes and zscores compared to the selected references
#'

rpm <- apply(as.matrix(reference[,choice$reference.choice]), 2, function(x) x / sum(x))

bg <- tibble(
    key = reference$key, 
    `mean reference` = rowMeans(rpm), 
    `sd reference` = apply(rpm, 1, sd))

amplicon <- left_join(block, bg, by = "key") %>%
    mutate(
        rpm = depth / sum(depth),
        zscore = (rpm - `mean reference`) / `sd reference`,
        fold = log2(rpm / `mean reference`))

matcher <- str_match(amplicon$key, "^(.+):([0-9]+)-([0-9]+)\\(([+,-])\\)")
amplicon <- amplicon %>%
    mutate(
        chromosome = matcher[,2],
        start = parse_integer(matcher[,3]),
        end = parse_integer(matcher[,4]),
        strand = matcher[,5])

#' write the amplicon statistics
write_tsv(amplicon, fn_amplicon)


# 
# tf <- with(amplicon, chromosome == "chr14" & start >= 105000000 & end <= 108000000)
# tx <- with(cnv_calls, chromosome == "14" & start >= 105000000 & end <= 108000000)
# 
# amplicon[tf,] %>%
#   ggplot(aes(x=start/1000000, y = fold)) +
#   geom_point(size = 1) +
#   geom_segment(aes(x = start/1000000, xend = end/1000000, y = 1, yend = 1), colour = "red", data = cnv_calls[tx, ]) +
#   theme_bw() +
#   xlab("Genome coordinate (Mb)")
# ggsave("figures/NA15510_chr14_CNV.png", width = 8, height = 3, dpi = 300)

