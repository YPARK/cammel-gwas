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

med.stat.tab <- bind_rows(lapply(med.files, read_tsv)) %>%
    separate(med.id, into = c('med.id', 'factor'), sep = '@') %>%
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
            mutate(pval.lodds = p.lodds.null)

    return(ret)
}

take.lfsr <- function(tab, var.min = 1e-4) {
    ret <- mutate(.data = tab, pip = 1/(1+exp(-lodds))) %>%
        mutate(pos.prob = pnorm(0, mean = theta, sd = sqrt(theta.var + var.min))) %>%
            mutate(neg.prob = 1 - pos.prob) %>%
                mutate(lfsr = 1 - pip * pmax(pos.prob, neg.prob)) %>%
                    select(-pos.prob, -neg.prob, -pip)
    return(ret)
}

gwas.names <- unique(med.stat.tab$gwas)

pval.tab <- bind_rows(lapply(gwas.names, calc.pval))

out.tab <- med.stat.tab %>% take.lfsr() %>%
    left_join(pval.tab)

write_tsv(out.tab, path = out.file)
