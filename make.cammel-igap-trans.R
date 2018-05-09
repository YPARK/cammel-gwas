#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-cammel.R')
library(readr)
library(dplyr)
library(zqtl)
library(methods)

if(length(argv) != 5) {
    q()
}

trans.eqtl.file <- argv[1]          # e.g., trans.eqtl.file <- 'trans-eqtl/igap_rosmap_gammax-4_eigen-2/chr19/1607_105.txt.gz'
geno.dir <- argv[2]                 # e.g., geno.dir <- 'ROSMAP_GENO'
gammax.input <- as.numeric(argv[3]) # e.g., gammax.input = 1e4
eig.tol <- as.numeric(argv[4])      # e.g., eig.tol = 1e-2
out.hdr <- argv[5]                  # e.g., out.hdr = 'temp.cammel'

################################################################
dir.create(dirname(out.hdr), recursive = TRUE)

trans.info <- (trans.eqtl.file %>% strsplit(split = '[/]'))[[1]]

data.info <- basename(trans.eqtl.file) %>%
    strsplit(split = '_') %>%
        .unlist()

gwas.ld.idx <- data.info[1] %>% as.integer()
lodds.cutoff <- log(0.9) - log(0.1)

trans.out.file <- out.hdr %&&% '.trans-mediation.txt.gz'

.files <- c(trans.out.file)
if(all(sapply(.files, file.exists))) {
    log.msg('all the output files exist: %s\n', paste(.files, collapse = ', '))
    q()
}

gwas.dir <- './gwas_stat/'
igap.gwas.file <- gwas.dir %&&% '/igap_' %&&% gwas.ld.idx %&&% '.txt.gz'

.files <- c(igap.gwas.file)
if(!all(sapply(.files, file.exists))) {
    log.msg('Missing input files: %s\n', paste(.files, collapse = ', '))
    q()
}

ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)
chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[gwas.ld.idx, 'chr']) %>%
    as.integer()
ld.lb.input <- ld.info.tab[gwas.ld.idx, 'start'] %>% as.integer()
ld.ub.input <- ld.info.tab[gwas.ld.idx, 'stop'] %>% as.integer()

################################################################

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(geno.dir %&&% '/chr' %&&% chr.input %&&% '.bed')) {
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
    write_tsv(data.frame(), path = trans.out.file)
    q()
}

plink.gwas <-  plink.matched$gwas
plink.eqtl <-  plink.matched$qtl

###########################
## Test 1 : G -> T -> AD ##
###########################

################################################################
## chr.snp rs snp.loc chr.gene med.id qtl.a1 qtl.a2 qtl.beta qtl.z
## i       c  i       i        c      c      c      d        d
trans.eqtl.tab <- read_tsv(trans.eqtl.file, col_types = 'iciicccdd')

qtl.tab <- trans.eqtl.tab %>%
    dplyr::rename(chr = chr.snp) %>%
        dplyr::select(-rs)

if(nrow(qtl.tab) == 0) {
    write_tsv(data.frame(), path = trans.out.file)
    log.msg('Empty QTL file\n')
    q()
}

igap.gwas.tab <- read_tsv(igap.gwas.file)

igap.data <- igap.gwas.tab %>%
    match.allele(plink.obj = plink.eqtl, qtl.tab = qtl.tab) %>%
        make.zqtl.data()
gc()

vb.opt <- list(pi.ub = -1, pi.lb = -5, tau = -5, do.hyper = TRUE, tol = 1e-8,
               gammax = gammax.input,
               vbiter = 3500, do.stdize = TRUE, eigen.tol = eig.tol,
               rate = 1e-2, nsample = 10, print.interv = 500,
               weight = FALSE, do.rescale = FALSE)

trans.out.tab <- igap.data %>%
    run.cammel(xx.gwas = plink.gwas$BED, xx.med = plink.eqtl$BED, opt = vb.opt) %>%
        get.effect.tab(z.data = igap.data,
                       gwas.tab = igap.gwas.tab,
                       qtl.tab = qtl.tab,
                       data.name = 'IGAP')

write_tsv(trans.out.tab, path = trans.out.file)

sig.med.id <- trans.out.tab %>%
    dplyr::filter(lodds > lodds.cutoff) %>%
        dplyr::select(med.id) %>% unique() %>% .unlist()

log.msg('Selected mediation IDs: %s\n', paste(sig.med.id, collapse = ', '))

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
