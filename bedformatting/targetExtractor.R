#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: targetExtractor.R
##
## Purpose of script: Extract targets from bed file
##
## Author: Niklas Reimer, M.Sc.
##
## Date Created: 2020-11-04
## Email: niklas.reimer@uksh.de
##
## ---------------------------

library(tidyr)
library(biomaRt)

args = commandArgs(trailingOnly=TRUE)

# read bed file
bed <- read.table(paste(args[1], ".bed", sep = ""))

# remove chr prefix if present
regions <- paste(gsub("chr", "", bed$V1), paste(bed$V2, bed$V3, sep="-"), sep=":")

# setup ensembl
m <- useMart('ensembl', dataset='hsapiens_gene_ensembl') # create a mart object
df <- getBM(mart=m, attributes=c('hgnc_symbol'), filters=c('chromosomal_region'), values=list(c(regions)))

# remove duplicate symbols
hgnc <- df$hgnc_symbol[!duplicated(df$hgnc_symbol)]

# write to file
sink(paste(args[1], "_Targets.txt", sep = ""))
cat(writeLines(hgnc, sep = "\n"))
sink()
