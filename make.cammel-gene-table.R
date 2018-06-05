#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(tidyr)
library(readr)
library(methods)

if(length(argv) < 3) {
    q()
}

med.dir <- argv[1] # e.g., med.dir = 'mediation.joint/ukbb/gammax_4/eigen_2'
med.ext <- argv[2] # e.g., med.ext = 'ad_mat'
out.file <- argv[3]

if(file.exists(out.file)) {
    log.msg('output file already exists')
    q()
}

all.genes <- read_tsv('all.genes.txt.gz',
                      col_names = c('chr', 'strand', 'med.id', 'hgnc', 'tss', 'tes'),
                      col_types = 'icccii',
                      skip = 1) %>%
                          na.omit()

library(biomaRt)

grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

coding.bm <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filter = 'biotype', values = 'protein_coding',
                   mart = grch37)

.list.files <- function(...) list.files(..., full.names = TRUE)

.unlist <- function(...) unlist(..., use.names = FALSE)

## Read the files

med.files <- .list.files(path = med.dir, pattern = med.ext %&&% '.gz')

ld.file <- 'ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed'

ld.tab <- read_tsv(ld.file) %>%
    rename(ld.lb = start, ld.ub = stop) %>%
        mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
            dplyr::filter(chr %in% as.character(1:22)) %>%
                mutate(chr = as.integer(chr))

med.stat.tab <- bind_rows(lapply(med.files, read_tsv)) %>%
    separate(med.id, into = c('med.id', 'factor', 'data'), sep = '@') %>%
        separate(med.id, into = c('med.id', 'remove'), sep = '[.]') %>%
            dplyr::select(-remove) %>%
                left_join(all.genes) %>%
                    dplyr::filter(hgnc %in% coding.bm$hgnc_symbol) %>%
                        group_by(med.id, factor, data) %>% slice(which.max(lodds)) %>%
                            dplyr::filter(!is.na(hgnc)) %>%
                                as.data.frame()

perm.stat.tab <- bind_rows(lapply(med.files, read_tsv)) %>%
    filter(substr(med.id, start = 2, stop = 5) == 'perm')

## local false sign rate
take.lfsr <- function(tab, var.min = 1e-8) {
    ret <- mutate(.data = tab, pip = 1/(1+exp(-lodds))) %>%
        mutate(pos.prob = pnorm(0, mean = theta, sd = sqrt(theta.var + var.min))) %>%
            mutate(neg.prob = 1 - pos.prob) %>%
                mutate(lfsr = 1 - pip * pmax(pos.prob, neg.prob)) %>%
                    dplyr::select(-pos.prob, -neg.prob, -pip) %>%
                        mutate(lfsr = signif(lfsr, 2))
    return(ret)
}

## Estimate empirical null distributions -- cohort by cohort
calc.pval <- function(tab, lodds.null) {
    null.mu <- mean(lodds.null)
    null.sd <- sd(lodds.null) + 1e-4
    pv.n <- pnorm((tab$lodds - null.mu) / null.sd, lower.tail = FALSE)

    cdf.0 <- ecdf(lodds.null)
    pv <- 1 - cdf.0(tab$lodds)
    tab %>% mutate(pval.emp = signif(pv, 2),
                   pval.nor = signif(pv.n, 2))
}

out.tab <- med.stat.tab %>% filter(num.genes.ld >= 20) %>%
    take.lfsr() %>%
        calc.pval(lodds.null = perm.stat.tab$lodds)

write_tsv(out.tab, path = out.file)
