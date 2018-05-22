#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(tidyr)
library(readr)
library(methods)

ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)

out.dir <- '/broad/hptmp/ypp/cammel/gwas_stat/'
system('mkdir -p ' %&&% out.dir)

## Divide GWAS statistics into independent chunks
write.gwas.chunk <- function(l, gwas.tab, gwas.name) {
    
    chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[l, 'chr']) %>%
        as.integer()
    ld.lb.input <- ld.info.tab[l, 'start'] %>% as.integer()
    ld.ub.input <- ld.info.tab[l, 'stop'] %>% as.integer()
   
    out.file <- out.dir %&&% '/' %&&% gwas.name %&&% '_' %&&% l %&&% '.txt.gz'

    if(!file.exists(out.file)) {
        gwas.tab %>% 
            dplyr::filter(chr == chr.input, snp.loc >= ld.lb.input, snp.loc <= ld.ub.input) %>%
                dplyr::mutate(gwas.z = gwas.beta / gwas.se) %>%
                    write_tsv(path = out.file)
        log.msg('Written %s', out.file)
        return(1)
    }
    return(0)
}

n.ld <- nrow(ld.info.tab)

## IGAP AD
igap.cols <- c('chr', 'snp.loc', 'rs', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p')
igap.gwas.file <- 'gwas/igap_ad.txt.gz'
igap.gwas.tab <- read_tsv(igap.gwas.file, col_names = igap.cols, skip = 1)

n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = igap.gwas.tab, gwas.name = 'igap')
log.msg('Wrote %d files\n', sum(n.written))
rm(igap.gwas.tab)
gc()

## other PGC gwas
pgc.cols <- c('chr', 'rs', 'snp.loc', 'gwas.a1', 'gwas.a2', 'info', 'gwas.beta', 'gwas.se', 'gwas.p')
pgc.tab <- read_tsv('gwas/adhd_jul2017.gz', col_names = pgc.cols, skip = 1) %>%
    na.omit() %>%
        mutate(chr = as.integer(chr), gwas.beta = log(gwas.beta))
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = pgc.tab, gwas.name = 'pgc_adhd')
log.msg('Wrote %d files\n', sum(n.written))
rm(pgc.tab)
gc()

pgc.cols <- c('chr', 'rs', 'snp.loc', 'gwas.a1', 'gwas.a2', 'info', 'gwas.beta', 'gwas.se', 'gwas.p')
pgc.tab <- read_tsv('gwas/ocd_aug2017.gz', col_names = pgc.cols, skip = 1) %>%
    na.omit() %>%
        mutate(chr = as.integer(chr), gwas.beta = log(gwas.beta))
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = pgc.tab, gwas.name = 'pgc_ocd')
log.msg('Wrote %d files\n', sum(n.written))
rm(pgc.tab)
gc()

pgc.cols <- c('chr', 'rs', 'gwas.a1', 'gwas.a2', 'snp.loc', 'info', 'gwas.beta', 'gwas.se', 'gwas.p', 'ngt')
pgc.tab <- read_tsv('gwas/ckqny.scz2snpres.gz', col_names = pgc.cols, skip = 1) %>%
    mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
        mutate(chr = as.integer(chr), gwas.beta = log(gwas.beta)) %>%
            na.omit()
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = pgc.tab, gwas.name = 'pgc_scz')
log.msg('Wrote %d files\n', sum(n.written))
rm(pgc.tab)
gc()

pgc.cols <- c('chr', 'rs', 'snp.loc', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p')
pgc.types <- 'icicc___ddd________'
pgc.tab <- read_tsv('gwas/pgc_mdd_tab.gz', col_names = pgc.cols, col_types = pgc.types, skip = 1) %>%
    na.omit() %>%
        mutate(chr = as.integer(chr), gwas.beta = log(gwas.beta))
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = pgc.tab, gwas.name = 'pgc_mdd')
log.msg('Wrote %d files\n', sum(n.written))
rm(pgc.tab)
gc()

pgc.cols <- c('chr', 'rs', 'snp.loc', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p')
pgc.types <- 'icicc___ddd_'
pgc.tab <- read_tsv('gwas/daner_pgc_asd_euro_all_25Mar2015.gz', col_names = pgc.cols, col_types = pgc.types, skip = 1) %>%
    na.omit() %>%
        mutate(chr = as.integer(chr), gwas.beta = log(gwas.beta))
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = pgc.tab, gwas.name = 'pgc_asd')
log.msg('Wrote %d files\n', sum(n.written))
rm(pgc.tab)
gc()

################################################################

read.ukbb <- function(...) {
    ret <- read_tsv(...) %>%
        separate('variant', c('chr', 'snp.loc', 'gwas.a1', 'gwas.a2')) %>%
            mutate(chr = as.integer(chr), snp.loc = as.integer(snp.loc)) %>%
                na.omit() %>%
                    rename(rs = rsid, gwas.beta = beta, gwas.se = se, gwas.p = pval) %>%
                        select(chr, snp.loc, rs, gwas.a1, gwas.a2, gwas.beta, gwas.se, gwas.p)
}

ukbb.paternal.gwas.file <- 'gwas/ukbb_paternal_ad_formatted.txt.gz'
ukbb.maternal.gwas.file <- 'gwas/ukbb_maternal_ad_formatted.txt.gz'

ukbb.pat.gwas.tab <- read.ukbb(ukbb.paternal.gwas.file)
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = ukbb.pat.gwas.tab, gwas.name = 'ukbb_ad_pat')
log.msg('Wrote %d files\n', sum(n.written))
rm(ukbb.pat.gwas.tab)
gc()

ukbb.mat.gwas.tab <- read.ukbb(ukbb.maternal.gwas.file)
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = ukbb.mat.gwas.tab, gwas.name = 'ukbb_ad_mat')
log.msg('Wrote %d files\n', sum(n.written))
rm(ukbb.mat.gwas.tab)
gc()

################################################################

ukbb.tab <- read.ukbb('gwas/ukbb_neuroticism_formatted.txt.gz')
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = ukbb.tab, gwas.name = 'ukbb_neuroticism')
log.msg('Wrote %d files\n', sum(n.written))
rm(ukbb.tab)
gc()

ukbb.tab <- read.ukbb('gwas/ukbb_moodswings_formatted.txt.gz')
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = ukbb.tab, gwas.name = 'ukbb_moodswings')
log.msg('Wrote %d files\n', sum(n.written))
rm(ukbb.tab)
gc()
