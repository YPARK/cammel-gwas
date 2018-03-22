#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(dplyr)
library(readr)
library(methods)
source('Util.R')

expr.file <- 'Mayo/MayoRNAseq_RNAseq_TCX_geneCounts.tsv'
covar.file <- 'Mayo/MayoRNAseq_RNAseq_TCX_covariates.csv'
fam.file <- 'Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed.fam'

all.genes <- read_tsv('all.genes.txt.gz',
                      col_names = c('chr', 'strand', 'ensembl_id', 'hgnc', 'tss', 'tes'),
                      col_types = 'icccii',
                      skip = 1) %>%
                          na.omit()

fam.cols <- c('FID', 'ID', 'paternal', 'maternal', 'sex', 'pheno')
fam.tab <- read_delim(fam.file, delim = ' ', col_names = fam.cols) %>%
    mutate(plink.idx = 1:n())

expr.data <- read_tsv(expr.file)

covar.tab <- read_csv(covar.file) %>%
    mutate(ID = gsub(ID, pattern = '_TCX', replacement = '')) %>%
        mutate(ID = as.integer(ID))

genes.info <- expr.data %>% select(ensembl_id) %>%
    mutate(idx = 1:n()) %>%
        left_join(all.genes) %>%
            na.omit()

sample.info <- tibble(ID = colnames(expr.data[, -1])) %>%
    mutate(ID = gsub(ID, pattern = '_TCX', replacement = '')) %>%
        mutate(ID = as.integer(ID)) %>%
            mutate(expr.idx = 2:(n() + 1)) %>%
                left_join(covar.tab) %>%
                    left_join(fam.tab %>% select(ID, plink.idx)) %>%
                        filter(!is.na(expr.idx), !is.na(plink.idx), !is.na(RIN))

Y <- expr.data[genes.info$idx, sample.info$expr.idx] %>% t()

genes.info <- genes.info %>% select(-idx)

.temp <- apply(Y < 1, 2, mean, na.rm = TRUE)
.rm <- which(is.na(.temp) | .temp > .5)

if(length(.rm) > 0) {
    Y <- Y[, -.rm, drop = FALSE]
    genes.info <- genes.info[ -.rm, , drop = FALSE]
}

Y <- Y %>% adjust.size.factor()

Y.out <- t(Y) %>% as.data.frame()
colnames(Y.out) <- sample.info$ID

dir.create('cis-eqtl/mayo_data', recursive = TRUE)

write_tsv(genes.info, path = 'cis-eqtl/mayo_data/genes.txt.gz')
write_tsv(sample.info, path = 'cis-eqtl/mayo_data/samples.txt.gz')
write_tsv(Y.out, path = 'cis-eqtl/mayo_data/expression.txt.gz')

