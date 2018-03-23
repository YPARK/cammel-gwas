#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(readr)
library(rtracklayer)

## select coding genes in GTF
gtf.file <- 'gencode.v19.genes.patched_contigs.gtf.gz'
out.file <- 'coding.genes.txt.gz'
out.file.2 <- 'all.genes.txt.gz'

gtf.tab <- readGFF(gtf.file, tags = c('gene_id', 'gene_name', 'transcript_name', 'gene_type'))

coding.genes <- gtf.tab %>% mutate(chr = seqid, ensg = gene_id) %>%
    filter(gene_type == 'protein_coding', type == 'transcript') %>%
        separate(ensg, into = c('ensg', 'remove'), sep = '[.]') %>%
            select(chr, start, end, strand, ensg, gene_name) %>%            
                arrange(chr, start)

all.genes <- gtf.tab %>% mutate(chr = seqid, ensg = gene_id) %>%
    filter(type == 'transcript') %>%
        separate(ensg, into = c('ensg', 'remove'), sep = '[.]') %>%
            group_by(chr, strand, ensg, gene_name) %>%            
                summarize(start = min(start), end = max(end)) %>%
                    arrange(chr, start)

write_tsv(coding.genes, path = out.file)

write_tsv(all.genes, path = out.file.2)
