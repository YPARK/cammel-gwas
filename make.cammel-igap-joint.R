#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

if(length(argv) != 4) {
    q()
}

ld.idx <- as.integer(argv[1])           # e.g., ld.idx = 133
gammax.input <- as.numeric(argv[2])     # e.g., gammax.input = 1e4
eig.tol <- as.numeric(argv[3])          # e.g., eig.tol = 1e-2
out.hdr <- argv[4]                      # e.g., out.hdr = 'temp.cammel'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

out.tab.file <- out.hdr %&&% '.mediation.gz'
out.null.file <- out.hdr %&&% '.null.gz'

.files <- c(out.tab.file, out.null.file)
if(all(sapply(.files, file.exists))) {
    log.msg('all the output files exist: %s\n', paste(.files, collapse = ', '))
    q()
}

################################################################
library(zqtl)
library(dplyr)
library(readr)
library(methods)
library(tidyr)

## Combine all eQTL effects
eqtl.data <- c('rosmap-mult', 'mayo-mult', 'geuvadis-mult', 'gtex-fqtl-v6')
eqtl.data.files <- 'cis-eqtl/' %&&% eqtl.data %&&% '/' %&&% ld.idx %&&% '_qtl.txt.gz'

################################################################
## Read QTL statistics and measure basic statistics
## chr rs snp.loc med.id qtl.a1 qtl.a2 qtl.beta qtl.z
## i   c  i       c      c      c      d        d
read.eqtl <- function(ii) {
    eqtl.cols <- 'icicccdd'
    ret <- read_tsv(eqtl.data.files[ii], col_types = eqtl.cols)
    if(nrow(ret) == 0) return(NULL)
    ret <- ret %>% mutate(data = eqtl.data[ii])
    return(ret)
}

eqtl.tab <- 1:length(eqtl.data) %>% lapply(FUN = read.eqtl) %>%
    bind_rows()

if(nrow(eqtl.tab) < 1) {
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = out.null.file)
    q()
}

eqtl.tab <- eqtl.tab %>%
    separate(col = med.id, into = c('med.id', 'factor'), sep = '@') %>%
        mutate(factor = if_else(is.na(factor), '0', factor)) %>%
            mutate(factor = as.integer(factor)) %>%
                separate(col = med.id, into = c('med.id', 'remove'), sep = '[.]') %>%
                    select(-remove)

eqtl.tab <- eqtl.tab %>%
    mutate(med.id = med.id %&&% '@' %&&% factor %&&% '@' %&&% data) %>%
        select(-factor, -data)

source('Util-cammel.R')

################################################################
ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)
chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[ld.idx, 'chr']) %>%
    as.integer()
ld.lb.input <- ld.info.tab[ld.idx, 'start'] %>% as.integer()
ld.ub.input <- ld.info.tab[ld.idx, 'stop'] %>% as.integer()

gwas.dir <- './gwas_stat/'
igap.gwas.file <- gwas.dir %&&% '/igap_' %&&% ld.idx %&&% '.txt.gz'

.files <- c(igap.gwas.file)
if(!all(sapply(.files, file.exists))) {
    log.msg('Missing input files: %s\n', paste(.files, collapse = ', '))
    q()
}

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

plink.gwas <- subset.plink('1KG_EUR/chr' %&&% chr.input,
                           chr.input, ld.lb.input, ld.ub.input, temp.dir)


if(nrow(eqtl.tab) == 0) {
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = out.null.file)
    log.msg('Empty QTL file\n')
    q()
}

igap.gwas.tab <- read_tsv(igap.gwas.file)

################################################################
igap.matched <- igap.gwas.tab %>%
    match.allele(plink.obj = plink.gwas, qtl.tab = eqtl.tab)

igap.data <- igap.matched %>%
    make.zqtl.data()
gc()

vb.opt <- list(pi.ub = -1, pi.lb = -5, tau = -5,
               do.hyper = TRUE, tol = 1e-8,
               gammax = gammax.input, nsingle = 150,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = eig.tol,
               rate = 1e-2, nsample = 10, print.interv = 500,
               weight = FALSE, do.rescale = FALSE)

if(!file.exists(out.tab.file)) {
    out.tab <- igap.data %>%
        run.cammel(xx.gwas = plink.gwas$BED, xx.med = plink.gwas$BED, opt = vb.opt) %>%
            get.effect.tab(igap.data, igap.gwas.tab, eqtl.tab, 'IGAP')
    
    out.tab <- out.tab %>% mutate(ld.idx = ld.idx)
    write_tsv(out.tab, path = out.tab.file)
}

if(!file.exists(out.null.file)) {
    null.tab <- igap.data %>%
        run.cammel.null(xx.gwas = plink.gwas$BED, xx.med = plink.gwas$BED, n.null = 1, opt = vb.opt) %>%
            mutate(gwas = 'IGAP')

    write_tsv(null.tab, path = out.null.file)
}

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
