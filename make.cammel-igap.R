#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

if(length(argv) != 6) {
    q()
}

ld.idx <- as.integer(argv[1])           # e.g., ld.idx = 133
qtl.file <- argv[2]                     # e.g., qtl.file = 'cis-eqtl/geuvadis/133_qtl.txt.gz'
geno.dir <- argv[3]                     # e.g., geno.dir = '1KG_EUR' # (eQTL genotype matrix)
gammax.input <- as.numeric(argv[4])     # e.g., gammax.input = 1e4, 1e3, 1e2
eig.tol <- as.numeric(argv[5])          # e.g., eig.tol = 1e-1
out.hdr <- argv[6]                      # e.g., out.hdr = 'temp.cammel'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

out.tab.file <- out.hdr %&&% '.mediation.gz'
out.null.file <- out.hdr %&&% '.null.gz'

.files <- c(out.tab.file, out.null.file)
if(all(sapply(.files, file.exists))) {
    log.msg('all the output files exist: %s\n', paste(.files, collapse = ', '))
    q()
}

if(!file.exists(qtl.file)) {
    log.msg('QTL file does not exist: %s\n', qtl.file)
    q()
}

################################################################
library(zqtl)
library(dplyr)
library(readr)
library(methods)

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

if(file.exists(geno.dir %&&% '/chr' %&&% chr.input %&&% '.bed')) {
    plink.eqtl <- subset.plink(geno.dir %&&% '/chr' %&&% chr.input,
                               chr.input, ld.lb.input, ld.ub.input, temp.dir)
} else if (file.exists(geno.dir %&&% '.bed')) {
    plink.eqtl <- subset.plink(geno.dir,
                               chr.input, ld.lb.input, ld.ub.input, temp.dir)
} else {
    plink.eqtl <- NULL
}

plink.gwas <- subset.plink('1KG_EUR/chr' %&&% chr.input,
                           chr.input, ld.lb.input, ld.ub.input, temp.dir)

## Read and match two PLINK filesets
plink.matched <- match.plink(plink.gwas, plink.eqtl)

if(is.null(plink.matched)) {
    log.msg('No common variants between eQTL and GWAS reference panels')
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = out.null.file)
    q()
}

plink.gwas <-  plink.matched$gwas
plink.eqtl <-  plink.matched$qtl

################################################################
## Read QTL statistics and measure basic statistics
## chr rs snp.loc med.id qtl.a1 qtl.a2 qtl.beta qtl.z
## i   c  i       c      c      c      d        d
clqtl.tab <- read_tsv(qtl.file, col_types = 'icicccdd')

if(nrow(qtl.tab) == 0) {
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = out.null.file)
    log.msg('Empty QTL file\n')
    q()
}

igap.gwas.tab <- read_tsv(igap.gwas.file)

################################################################
igap.matched <- igap.gwas.tab %>%
    match.allele(plink.obj = plink.eqtl, qtl.tab = qtl.tab)

igap.data <- igap.matched %>%
    make.zqtl.data()
gc()

vb.opt <- list(pi.ub = -1, pi.lb = -5, tau = -5, do.hyper = TRUE, tol = 1e-8,
               gammax = gammax.input,
               vbiter = 3500, do.stdize = TRUE, eigen.tol = eig.tol,
               rate = 1e-2, nsample = 10, print.interv = 500,
               weight = FALSE, do.rescale = FALSE)

out.tab <- igap.data %>%
    run.cammel(xx.gwas = plink.gwas$BED, xx.med = plink.eqtl$BED, opt = vb.opt) %>%
        get.effect.tab(z.data = igap.data, gwas.tab = igap.gwas.tab, qtl.tab = qtl.tab, data.name = 'IGAP')

write_tsv(out.tab, path = out.tab.file)

null.tab <- igap.data %>%
    run.cammel.null(xx.gwas = plink.gwas$BED, xx.med = plink.eqtl$BED, n.null = 5, opt = vb.opt) %>%
        mutate(gwas = 'IGAP')

write_tsv(null.tab, path = out.null.file)

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
