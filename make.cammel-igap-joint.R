#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-cammel.R')

if(length(argv) < 5) {
    q()
}

ld.idx <- as.integer(argv[1])           # e.g., ld.idx = 1
gammax.input <- as.numeric(argv[2])     # e.g., gammax.input = 1e4
eig.tol <- as.numeric(argv[3])          # e.g., eig.tol = 1e-2
eqtl.data <- argv[4]                    # e.g., eqtl.data = 'rosmap-mult:mayo-mult:geuvadis-mult'
out.hdr <- argv[5]                      # e.g., out.hdr = 'temp.cammel'

is.null.data <- FALSE
if(length(argv) > 5) {
    is.null.data <- as.logical(argv[6])
}

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

out.tab.file <- out.hdr %&&% '.igap_ad.gz'
z.tab.file <- out.hdr %&&% '.zscores.gz'

.files <- c(out.tab.file, z.tab.file)
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
eqtl.data <- strsplit(eqtl.data, split = ':') %>% (function(x) x[[1]])
eqtl.data.files <- 'cis-eqtl/' %&&% eqtl.data %&&% '/' %&&% ld.idx %&&% '_qtl.txt.gz'
eqtl.tab <- read.multivar.eqtl(eqtl.data, eqtl.data.files)

if(nrow(eqtl.tab) == 0) {
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = z.tab.file)
    log.msg('Empty QTL file\n')
    q()
}

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

igap.gwas.tab <- read.gwas(igap.gwas.file)

################################################################
igap.matched <- igap.gwas.tab %>%
    match.allele(plink.obj = plink.gwas, qtl.tab = eqtl.tab)

igap.data <- igap.matched %>%
    make.zqtl.data(n.permuted = 50, is.null.data = is.null.data)
gc()

if(is.null(igap.data)) {
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = z.tab.file)
    log.msg('Empty data\n')
    q()
}

if(!file.exists(out.tab.file)) {

    vb.opt <- list(pi = -1, tau = -5, do.hyper = FALSE, vbiter = 5000, tol = 1e-8,
                   gammax = gammax.input, nsingle = 100,
                   do.stdize = TRUE, eigen.tol = eig.tol,
                   rate = 1e-2, decay = -1e-2, nsample = 11, print.interv = 500,
                   weight = FALSE, do.rescale = TRUE, jitter = 0.1,
                   multivar.mediator = TRUE)

    z.out <- igap.data %>%
        run.cammel(xx.gwas = plink.gwas$BED, xx.med = plink.gwas$BED, opt = vb.opt)

    var.tab <- get.var.tab(z.out$var.decomp, igap.data$mediators) %>%
        select(med.id, var.mediated, var.direct.tot)

    summary.tab <- igap.matched %>% group_by(med.id) %>% slice(which.max(abs(qtl.z))) %>%
        select(med.id, gwas.p, gwas.z)

    out.tab <- melt.effect(z.out$param.mediated, igap.data$mediators, 'IGAP') %>%
        rename(med.id = Var1, gwas = Var2) %>%
            left_join(var.tab) %>%
                left_join(summary.tab)

    out.tab <- out.tab %>%
        mutate(ld.idx = ld.idx,
               gwas.p.ld = min(igap.matched$gwas.p),
               num.genes.ld = nrow(summary.tab))

    if(sum(out.tab$lodds > 0) > 0) {
        zscore.tab <- separate.zscore(z.out, plink.gwas$BIM %r% igap.data$x.pos)
        write_tsv(zscore.tab, path = z.tab.file)
    } else {
        write_tsv(data.frame(), path = z.tab.file)
    }
    write_tsv(out.tab, path = out.tab.file)
}

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
