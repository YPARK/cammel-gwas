#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

if(length(argv) != 3) {
    q()
}

med.file <- argv[1]           # e.g., med.file = 'mediation/gene_igap_rosmap_gammax-4_eigen-2.txt.gz'
cutoff <- as.numeric(argv[2]) # e.g., cutoff = 0.01
out.file <- argv[3]           # e.g., ''

library(dplyr)
library(readr)
library(tidyr)
source('Util.R')

ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file, col_names = c('chr', 'ld.lb', 'ld.ub'), skip = 1) %>%
    mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
        mutate(chr = as.integer(chr)) %>%
            mutate(ld = chr %&&% ':' %&&% ld.lb %&&% ':' %&&% ld.ub) %>%
                mutate(idx = 1:n()) %>%
                    select(ld, idx)

dir.create(dirname(out.file), recursive = TRUE)

################################################################
med.tab <- read_tsv(med.file)

hgnc.ld.tab <- med.tab %>% select(chr, ld.lb, ld.ub, hgnc) %>% unique() %>%
    mutate(ld = chr %&&% ':' %&&% ld.lb %&&% ':' %&&% ld.ub) %>%
        select(hgnc, ld)

net.file <- 'pathways/BIOGRID-ORGANISM-Homo_sapiens-3.4.159.tab2.txt.gz'
net.tab <- read_tsv(net.file) %>%
    dplyr::filter(`Experimental System Type` == 'physical') %c%
    c(8, 9) %>% unique()
colnames(net.tab) <- c('gene1', 'gene2')

################################################################
## 1. Take significant genes.  These genes will be seed genes.  And
## construct subnetworks from the seed, looking for immediate
## neighbors.
sig.hgnc <- med.tab %>% dplyr::filter(lfsr < cutoff) %>%
    dplyr::select(hgnc) %>% unique() %>% .unlist()

net.genes <- net.tab %>% dplyr::filter(gene1 %in% sig.hgnc | gene2 %in% sig.hgnc) %>%
    dplyr::select(gene1, gene2) %>% .unlist() %>% unique()

net.genes <- c(net.genes, sig.hgnc) %>% unique() %>% sort()

net.edges <- net.tab %>% dplyr::filter(gene1 %in% net.genes, gene2 %in% net.genes) %>%
    dplyr::filter(gene1 != gene2) %>%
    unique()

net.flip <- net.edges %>% mutate(temp = gene2) %>%
    mutate(gene2 = gene1) %>%
        mutate(gene1 = temp) %>%
            select(-temp)

net.edges <- bind_rows(net.flip, net.edges) %>%
    filter(gene1 %in% sig.hgnc)

################################################################
out <- net.edges %>%
    left_join(hgnc.ld.tab %>% rename(gene1 = hgnc), by = 'gene1') %>%
        left_join(hgnc.ld.tab %>% rename(gene2 = hgnc), by = 'gene2', suffix = c('.src','.tgt')) %>%
            na.omit() %>%
                select(ld.src, ld.tgt) %>%
                    filter(ld.src != ld.tgt) %>%
                        unique()

out <-
    out %>% left_join(ld.info.tab %>% rename(ld.src = ld, idx.src = idx)) %>%
        left_join(ld.info.tab %>% rename(ld.tgt = ld, idx.tgt = idx)) %>%
            arrange(idx.src)

write_tsv(out, path = out.file, col_names = FALSE)
