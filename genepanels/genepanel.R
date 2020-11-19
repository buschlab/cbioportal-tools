#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: bed2cBioPortalGenePanel.R
##
## Converts a bed file to a gene panel for cBioPortal.
##
## Author: Niklas Reimer, M.Sc.
##
## Date Created: 2020-10-30
## Email: niklas.reimer@uksh.de
##
## ---------------------------

require(tidyr)
require(biomaRt)

args = commandArgs(trailingOnly=TRUE)

# read genes provided by cBioPortal
cbiogenes <- read.table("cbioportal_genes.txt", header = TRUE, sep = "\t", quote="\"")
# read provided metadata for panel
meta <- read.table(paste0(args[1], "/meta.csv"), header = TRUE, sep=";")
# read provided bed file
bed <- read.table(paste0(args[1], "/target.bed"))

# remove chr prefix if present
regions <- paste(gsub("chr", "", bed$V1), paste(bed$V2, bed$V3, sep="-"), sep=":")

# setup ensembl
m <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
df <- getBM(mart=m, attributes=c('entrezgene_id'), filters=c('chromosomal_region'), values=list(c(regions)))

# remove duplicate ids
entrez <- df$entrezgene_id[!duplicated(df$entrezgene_id)]

#get genes using entrez id from cbioportal dataset
genes <- cbiogenes$HUGO_GENE_SYMBOL[!is.na(match(cbiogenes$ENTREZ_GENE_ID, entrez))]

#sort genes
genes <- sort(genes)

# write to file
sink(paste0(args[1], "/panel.txt"))
cat(paste("stable_id", meta$stable_id, sep = ": "))
cat("\n")
cat(paste("description", meta$description, sep = ": "))
cat("\n")

cat("gene_list: ")
cat(writeLines(as.character(genes), sep="\t"))
sink()
