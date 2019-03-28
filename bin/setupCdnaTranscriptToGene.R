#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(rtracklayer))

cl <- commandArgs(trailingOnly = TRUE)

cdna_file <- cl[1]
gtf_file <- cl[2]
transcript_to_gene_file <- cl[3]
filtered_cdna_file <- cl[4]

# Derive a transcript-to-gene mapping from the GTF

annotation <- elementMetadata(import(gtf_file))
annotation <- subset(annotation, ! (is.na(transcript_id) | is.na(gene_id)))
mapping <- cbind(paste(annotation$transcript_id, annotation$transcript_version, sep='.'), annotation$gene_id)

# Read the cDNA

cdna <- readDNAStringSet(cdna_file)
cdna_transcript_names <- unlist(lapply(names(cdna), function(x) unlist(strsplit(x, ' '))[1]  ))

# Filter out cDNAs without matching transcript entries in the GTF

cdna <- cdna[which(cdna_transcript_names %in% mapping[,1])]

# Write transcript-to-gene mapping and filtered cDNA file

write.table(unique(mapping), file= transcript_to_gene_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
writeXStringSet(x = cdna, filepath = filtered_cdna_file, compress = 'gzip')
