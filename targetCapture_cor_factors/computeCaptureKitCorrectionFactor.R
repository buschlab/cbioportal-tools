#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: computeCaptureKitCorrectionFactor.R
##
## Purpose of script: Compure correction factors for triplets between WGS and WES pipelines.
## For details you can have a look at the documentation of the R package YAPSA:
## https://www.bioconductor.org/packages/release/bioc/vignettes/YAPSA/inst/doc/vignette_exomes.html#introduction
##
## Author: Niklas Reimer, M.Sc.
##
## Date Created: 2020-11-04
## Email: niklas.reimer@uksh.de
##
## ---------------------------

require(bedr)
require(plyr)
require(stringr)

args = commandArgs(trailingOnly=TRUE)

# read bed file
bed <- read.table(args[1])
bed <- transform(bed, V1 = as.character(V1), V2 = as.numeric(V2), V3 = as.numeric(V3))
bed <- bedr.sort.region(bed)

# extract regions from fasta
fasta <- get.fasta(bed, fasta = args[2])

# load targetCapture_cor_factors
load(args[3])

# read triplets with every possible shift
tripletlist1 <- strsplit(as.character(fasta$sequence), split = "(?<=.{3})", perl = TRUE)
tripletlist2 <- strsplit(as.character(substring(fasta$sequence, 2)), split = "(?<=.{3})", perl = TRUE)
tripletlist3 <- strsplit(as.character(substring(fasta$sequence, 3)), split = "(?<=.{3})", perl = TRUE)
triplets1 <- unlist(tripletlist1)
triplets2 <- unlist(tripletlist2)
triplets3 <- unlist(tripletlist3)
triplets <- c(triplets1, triplets2, triplets3)

# filter by sequence
counts <- count(toupper(triplets))
# remove sequences with length lower than 3
countsFiltered <- counts[str_length(counts$x) == 3,]

# remove ambigious triplets
countsFiltered <- countsFiltered[str_detect(countsFiltered$x, "N", negate = TRUE),]

# calculate correction factors
abs_cor <- targetCapture_cor_factors["hs37d5"]$hs37d5$abs_cor / countsFiltered$freq
rel_cor <- abs_cor / (sum(targetCapture_cor_factors["hs37d5"]$hs37d5$abs_cor) / sum(countsFiltered$freq))

# save computated data
targetCapture_cor_factors[[args[4]]] <- list("rel_cor" = rel_cor, "abs_cor" = abs_cor)
save(targetCapture_cor_factors, file = file.path(args[5]))
