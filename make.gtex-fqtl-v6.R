#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

if(length(argv) < 4) {
    q()
}

chr.input <- as.integer(argv[1])   # chr.input = 1
ld.lb.input <- as.integer(argv[2]) # ld.lb.input = 11777841
ld.ub.input <- as.integer(argv[3]) # ld.ub.input = 12779466
out.file <- argv[4]                # out.file = 'temp.txt.gz'

if(file.exists(out.file)) {
    log.msg('File exists: %s\n', out.file)
    q()
}

geno.dir <- 'GTEX_GENO'
gtex.stat.file <- 'GTEX_v6_FQTL/chr' %&&% chr.input %&&% '/50/combined.txt.gz'

cis.dist <- 1e6

library(fqtl)
library(dplyr)
library(tidyr)
library(readr)
library(methods)

################################################################
## Read multivariate effect sizes
gtex.cols <- c('ensg', 'chr', 'tss',
               'tis.idx', 'tis.names', 'tis.theta', 'tis.se', 'tis.lodds',
               'snp.names', 'snp.theta', 'snp.se', 'snp.lodds',
               'factor', 'pip.cutoff')

.split.bar <- function(s) strsplit(s, split = '[|]')

gtex.stat <- read_tsv(gtex.stat.file, col_names = gtex.cols) %>%
    filter(pip.cutoff == .5) %>%
        filter(tss > (ld.lb.input - cis.dist), tss < (ld.ub.input + cis.dist))

if(nrow(gtex.stat) < 1) {
    log.msg('No valid gene\n')
    write_tsv(data.frame(), path = out.file)
    q()
}

gtex.stat <- gtex.stat %>%
    mutate(med.id = ensg %&&% '@' %&&% factor) %>%
        select(med.id, snp.names, snp.theta, snp.se, snp.lodds)  %>%
            unnest(snp = .split.bar(snp.names), theta = .split.bar(snp.theta), theta.se = .split.bar(snp.se), lodds = .split.bar(snp.lodds)) %>%
                select(-snp.names, -snp.theta)

snp.cols <- c('chr', 'snp.loc', 'a1', 'a2', 'remove')

gtex.stat <-
    gtex.stat %>%
        mutate(rs = snp) %>%
            separate(snp, snp.cols, sep = '[_]') %>%
                select(-remove) %>%
                    mutate(snp.loc = as.integer(snp.loc),
                           theta = as.numeric(theta),
                           theta.se = as.numeric(theta.se),
                           lodds = as.numeric(lodds),
                           chr = as.integer(chr))

## compare two plink bim files
read.bim <- function(plink.hdr, plink.lb, plink.ub) {
    bim <- read_tsv(plink.hdr %&&% '.bim',
                    col_names = c('chr', 'rs', 'missing', 'snp.loc', 'a1', 'a2'),
                    col_types = 'iciicc')

    bim <- bim %>% filter(snp.loc >= plink.lb, snp.loc <= plink.ub)
    return(bim)
}

match.bim <- function(gwas.bim, qtl.bim) {

    gwas.bim <- gwas.bim %>%
        mutate(gwas.x.pos = 1:n()) %>%
            rename(gwas.plink.a1 = a1,
                   gwas.plink.a2 = a2) %>%
                       select(-missing)

    qtl.bim <- qtl.bim %>%
        mutate(qtl.x.pos = 1:n()) %>%
            rename(qtl.plink.a1 = a1,
                   qtl.plink.a2 = a2,
                   qtl.rs = rs) %>%
                       select(-missing)

    bim.matched <- gwas.bim %>%
        left_join(qtl.bim) %>%
            na.omit()

    return(bim.matched)
}

gtex.bim <- (geno.dir %&&% '/chr' %&&% chr.input) %>%
    read.bim(plink.lb = ld.lb.input, plink.ub = ld.ub.input)
gwas.bim <- ('1KG_EUR/chr' %&&% chr.input) %>%
    read.bim(plink.lb = ld.lb.input, plink.ub = ld.ub.input)
matched.bim <- match.bim(gwas.bim, gtex.bim) %>%
    filter(gwas.plink.a1 == qtl.plink.a1, gwas.plink.a2 == qtl.plink.a2)

if(nrow(matched.bim) == 0) {
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)    
    q()
}

################################################################
gtex.snp.stat <- gtex.stat %>%
    select(med.id, snp.loc, a1, a2, theta, theta.se, lodds)

qtl.out <- gtex.snp.stat %>% left_join(matched.bim) %>% na.omit() %>%
    select(-qtl.rs) %>%
        rename(qtl.beta = theta, qtl.se = theta.se, qtl.lodds = lodds)

qtl.out <- qtl.out %>%
    mutate(qtl.beta = if_else(a2 != gwas.plink.a1, -qtl.beta, qtl.beta))    

qtl.out <- qtl.out %>%
    rename(qtl.a1 = gwas.plink.a1, qtl.a2 = gwas.plink.a2) %>%
        select(chr, rs, snp.loc, med.id, qtl.a1, qtl.a2, qtl.beta, qtl.se, qtl.lodds)

write_tsv(x = qtl.out, path = out.file)
system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
