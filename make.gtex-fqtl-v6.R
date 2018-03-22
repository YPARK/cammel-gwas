#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

if(length(argv) < 4) {
    q()
}

chr.input <- as.integer(argv[1])   # chr.input = 1
ld.lb.input <- as.integer(argv[2]) # ld.lb.input = 9365199
ld.ub.input <- as.integer(argv[3]) # ld.ub.input = 10806984
out.file <- argv[4]

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
        select(med.id, snp.names, snp.theta)  %>%
            unnest(snp = .split.bar(snp.names), theta = .split.bar(snp.theta)) %>%
                select(-snp.names, -snp.theta)

snp.cols <- c('chr', 'snp.loc', 'plink.a1', 'plink.a2', 'remove')

gtex.stat <-
    gtex.stat %>%
        mutate(rs = snp) %>%
            separate(snp, snp.cols, sep = '[_]') %>%
                select(-remove) %>%
                    mutate(snp.loc = as.integer(snp.loc),
                           theta = as.numeric(theta),
                           chr = as.integer(chr))

## Read genotype matrix and expand it to univariate statistics

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel.gtex-fqtl/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel.gtex-fqtl/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

plink.hdr <- geno.dir %&&% '/chr' %&&% chr.input

plink.error <- function(e) {
    log.msg('No QTL here!\n')
    return(NULL)
}

plink <- tryCatch(subset.plink(plink.hdr, chr.input, ld.lb.input, ld.ub.input, temp.dir), error = plink.error)

if(is.null(plink)){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)    
    q()
}

x.bim <- plink$BIM %>%
    mutate(x.pos = 1:n()) %>%
        select(-missing)

gtex.stat.bim <- gtex.stat %>%
    left_join(x.bim) %>%
        na.omit()

if(nrow(gtex.stat.bim) < 1) {
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)    
    q()
}

X <- plink$BED %>% scale() %>% rm.na.zero()

x.bim <- plink$BIM %>%
    mutate(x.col = 1:n())

multi.to.uni <- function(tab) {
    eta <- (X %c% tab$x.pos) %*% matrix(as.numeric(tab$theta), ncol = 1)
    ret <- calc.qtl.stat(X, eta) %>%
        left_join(x.bim) %>%
            rename(qtl.a1 = plink.a1, qtl.a2 = plink.a2) %>%
                mutate(qtl.beta = signif(beta, 4), qtl.z = signif(beta/se, 4))
}

qtl.out <- gtex.stat.bim %>% group_by(med.id) %>%
    do(stat = multi.to.uni(.)) %>%
        unnest() %>%
            select(chr, rs, snp.loc, med.id, dplyr::starts_with('qtl'))

write_tsv(x = qtl.out, path = out.file)
system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
