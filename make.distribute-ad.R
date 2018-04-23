#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(readr)
library(methods)

ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)

out.dir <- '/broad/hptmp/ypp/cammel/gwas_stat/'
system('mkdir -p ' %&&% out.dir)

igap.gwas.file <- 'gwas/igap_ad.txt.gz'

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

igap.cols <- c('chr', 'snp.loc', 'rs', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p')
igap.gwas.tab <- read_tsv(igap.gwas.file, col_names = igap.cols, skip = 1)

n.ld <- nrow(ld.info.tab)
n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = igap.gwas.tab, gwas.name = 'igap')
log.msg('Wrote %d files\n', sum(n.written))
rm(igap.gwas.tab)
gc()

## ukbb.paternal.gwas.file <- 'gwas/ukbb_paternal_ad_formatted.txt.gz'
## ukbb.maternal.gwas.file <- 'gwas/ukbb_maternal_ad_formatted.txt.gz'

## ukbb.cols <- c('chr', 'snp.loc', 'gwas.a1', 'gwas.a2', 'rs', 'sample.size', 'gwas.beta',
##                'gwas.se', 'gwas.p')

## ukbb.pat.gwas.tab <- read_tsv(ukbb.paternal.gwas.file, col_names = ukbb.cols)
## n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = ukbb.pat.gwas.tab, gwas.name = 'ukbb_ad_pat')
## log.msg('Wrote %d files\n', sum(n.written))
## rm(ukbb.pat.gwas.tab)
## gc()

## ukbb.mat.gwas.tab <- read_tsv(ukbb.maternal.gwas.file, col_names = ukbb.cols)
## n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = ukbb.mat.gwas.tab, gwas.name = 'ukbb_ad_mat')
## log.msg('Wrote %d files\n', sum(n.written))
## rm(ukbb.mat.gwas.tab)
## gc()

