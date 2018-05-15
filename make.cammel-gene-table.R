#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(tidyr)
library(readr)
library(methods)

if(length(argv) < 2) {
    q()
}

med.dir <- argv[1] # e.g., med.dir = 'mediation/ad/rosmap/gammax_3/eigen_1'
out.file <- argv[2]

if(file.exists(out.file)) {
    log.msg('output file already exists')
    q()
}

all.genes <- read_tsv('all.genes.txt.gz',
                      col_names = c('chr', 'strand', 'med.id', 'hgnc', 'tss', 'tes'),
                      col_types = 'icccii',
                      skip = 1) %>%
                          na.omit()

.list.files <- function(...) list.files(..., full.names = TRUE)

.unlist <- function(...) unlist(..., use.names = FALSE)

## Read the files

med.files <- .list.files(path = med.dir, pattern = 'mediation.gz')

null.files <- .list.files(path = med.dir, pattern = 'null.gz')

null.stat.tab <- bind_rows(lapply(null.files, read_tsv))

ld.file <- 'ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed'
ld.tab <- read_tsv(ld.file) %>%
    rename(ld.lb = start, ld.ub = stop) %>%
        mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
            filter(chr %in% as.character(1:22)) %>%
                mutate(chr = as.integer(chr))

.read.med <- function(.file) {
    ld.idx <- basename(.file) %>%
        gsub(pattern = '.mediation.gz', replacement = '') %>%
            as.integer()

    ret <- read_tsv(.file)

    if(!'ld.lb' %in% colnames(ret)) {
        ld.info <- ld.tab[ld.idx, ] %>% unlist()
        ret <- ret %>%
            mutate(chr = ld.info[1], ld.lb = ld.info[2], ld.ub = ld.info[3])
    }
    return(ret)
}

med.stat.tab <- bind_rows(lapply(med.files, .read.med)) %>%
    separate(med.id, into = c('med.id', 'factor', 'data'), sep = '@') %>%
        separate(med.id, into = c('med.id', 'remove'), sep = '[.]') %>%
            select(-remove) %>%
                left_join(all.genes) %>%
                    filter(!is.na(hgnc))

if(all(is.na(med.stat.tab$factor))) {
    med.stat.tab <- med.stat.tab %>% mutate(factor = '1')
}

## Estimate empirical null distributions -- cohort by cohort

calc.pval <- function(.gwas) {

    .null <- null.stat.tab %>% filter(gwas == .gwas) %>%
        select(lodds) %>%
        .unlist()

    null.cdf <- ecdf(1/(1+exp(-.null)))

    .pval.tab <- med.stat.tab %>% filter(gwas == .gwas) %>%
        select(chr, ld.lb, ld.ub, med.id, factor, gwas, lodds) %>%
            mutate(pip = 1/(1+exp(-lodds)))

    p.lodds.null <- 1 - null.cdf(.pval.tab$pip)

    ret <- .pval.tab %>%
        select(chr, ld.lb, ld.ub, med.id, factor, gwas) %>%
            mutate(pval.lodds = p.lodds.null) %>%
                mutate(pval.lodds = signif(pval.lodds, 2))

    return(ret)
}

take.lfsr <- function(tab, var.min = 1e-4) {
    ret <- mutate(.data = tab, pip = 1/(1+exp(-lodds))) %>%
        mutate(pos.prob = pnorm(0, mean = theta, sd = sqrt(theta.var + var.min))) %>%
            mutate(neg.prob = 1 - pos.prob) %>%
                mutate(lfsr = 1 - pip * pmax(pos.prob, neg.prob)) %>%                    
                    select(-pos.prob, -neg.prob, -pip) %>%
                        mutate(lfsr = signif(lfsr, 2))
    return(ret)
}

gwas.names <- unique(med.stat.tab$gwas)

pval.tab <- bind_rows(lapply(gwas.names, calc.pval))

out.tab <- med.stat.tab %>% take.lfsr() %>%
    left_join(pval.tab)

write_tsv(out.tab, path = out.file)
