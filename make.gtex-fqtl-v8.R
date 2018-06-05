#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

## 1. take FQTL results
## 2. liftover hg38 to hg19

if(length(argv) != 5) {
    q()
}

chr.input <- as.integer(argv[1])   # chr.input = 1
ld.lb.input <- as.integer(argv[2]) # ld.lb.input = 11777841
ld.ub.input <- as.integer(argv[3]) # ld.ub.input = 12779466
gtex.stat.file <- argv[4]          # e.g., 'GTEX_v8_FQTL/fqtl_' %&&% chr.input %&&% '.txt.gz'
out.file <- argv[5]                # out.file = 'temp.txt.gz'

options(stringsAsFactors = FALSE)
source('Util.R')

if(file.exists(out.file)) {
    log.msg('File exists: %s\n', out.file)
    q()
}

################################################################

geno.dir <- 'GTEX_GENO_v8'
cis.dist <- 1e6

################################################################

library(fqtl)
library(dplyr)
library(tidyr)
library(readr)
library(methods)

################################################################
## Read multivariate effect sizes
.split.bar <- function(s) strsplit(s, split = '[|]')

## hg19 coding genes
coding.genes.hg19 <- read_tsv('coding.genes.txt.gz') %>%
    rename(hgnc = gene_name, lb.hg19 = start, ub.hg19 = end) %>%
        na.omit()

## hg38 gtex stat
gtex.hg38 <- read_tsv(gtex.stat.file) %>%
    separate(ensg, c('ensg', 'remove')) %>%
        select(-remove) %>%
            mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
                mutate(chr = as.integer(chr)) %>%
                    left_join(coding.genes.hg19) %>%
                        na.omit()

snp.cols <- c('chr', 'snp.loc', 'a1', 'a2')

gtex.stat.hg38 <- gtex.hg38 %>%
    filter((lb.hg19 > (ld.lb.input - cis.dist) & lb.hg19 < (ld.ub.input + cis.dist)) | (ub.hg19 > (ld.lb.input - cis.dist) & ub.hg19 < (ld.ub.input + cis.dist)))

if(nrow(gtex.stat.hg38) < 1) {
    log.msg('No valid gene\n')
    write_tsv(data.frame(), path = out.file)
    q()
}

gtex.stat.hg38 <- gtex.stat.hg38 %>%
    mutate(med.id = ensg %&&% '@' %&&% factor) %>%
        select(med.id, snp, snp.theta, snp.sd, snp.lodds) %>%
            unnest(snp = .split.bar(snp),
                   snp.theta = .split.bar(snp.theta),
                   snp.sd = .split.bar(snp.sd),
                   snp.lodds = .split.bar(snp.lodds)) %>%
                       mutate(rs = snp) %>%
                           separate(snp, snp.cols, sep = ':')

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel-gtex/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel-gtex/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

gtex.stat <- gtex.stat.hg38 %>%
    lift.over.hg38.hg19(temp.hdr = temp.dir %&&% '/liftover')

gtex.bim <- (geno.dir %&&% '/chr' %&&% chr.input) %>%
    read.bim(plink.lb = ld.lb.input, plink.ub = ld.ub.input) %>%
        mutate(chr = 'chr' %&&% chr) %>%
            lift.over.hg38.hg19(temp.hdr = temp.dir %&&% '/liftover') %>%
                mutate(chr = as.integer(gsub(chr, pattern = 'chr', replacement = '')))

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
    mutate(chr = as.integer(gsub(chr, pattern = 'chr', replacement = ''))) %>%
        select(med.id, snp.loc, a1, a2, snp.theta, snp.sd, snp.lodds)

qtl.out <- gtex.snp.stat %>% left_join(matched.bim) %>% na.omit() %>%
    rename(qtl.beta = snp.theta, qtl.se = snp.sd, qtl.lodds = snp.lodds)

qtl.out <- qtl.out %>%
    mutate(qtl.beta = if_else(a2 != gwas.plink.a1, -as.numeric(qtl.beta), as.numeric(qtl.beta)))

qtl.out <- qtl.out %>%
    rename(qtl.a1 = gwas.plink.a1, qtl.a2 = gwas.plink.a2) %>%
        select(chr, rs, snp.loc, med.id, qtl.a1, qtl.a2, qtl.beta, qtl.se, qtl.lodds)

lodds.cutoff <- 0 # at least one pip > .5

valid.med <- qtl.out %>%
    group_by(med.id) %>%
        slice(which.max(qtl.lodds)) %>%
            filter(qtl.lodds > lodds.cutoff) %>%
                select(med.id) %>%
                    .unlist()

qtl.out <- qtl.out %>%
    filter(med.id %in% valid.med)

if(nrow(qtl.out) < 1){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

write_tsv(x = qtl.out, path = out.file)
system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
