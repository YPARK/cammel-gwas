#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-cammel.R')

if(length(argv) != 4) {
    q()
}

ld.idx <- as.integer(argv[1])           # e.g., ld.idx = 1
gammax.input <- as.numeric(argv[2])     # e.g., gammax.input = 1e4
eig.tol <- as.numeric(argv[3])          # e.g., eig.tol = 1e-2
out.hdr <- argv[4]                      # e.g., out.hdr = 'temp.cammel'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

################################################################
library(zqtl)
library(dplyr)
library(readr)
library(methods)
library(tidyr)

## Combine all eQTL effects
eqtl.data <- c('rosmap-mult', 'mayo-mult', 'geuvadis-mult', 'gtex-fqtl-v6')
eqtl.data.files <- 'cis-eqtl/' %&&% eqtl.data %&&% '/' %&&% ld.idx %&&% '_qtl.txt.gz'
eqtl.tab <- read.multivar.eqtl(eqtl.data, eqtl.data.files)
log.msg('Read %d rows eQTL tab\n', nrow(eqtl.tab))

################################################################
ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)
chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[ld.idx, 'chr']) %>%
    as.integer()
ld.lb.input <- ld.info.tab[ld.idx, 'start'] %>% as.integer()
ld.ub.input <- ld.info.tab[ld.idx, 'stop'] %>% as.integer()

gwas.dir <- './gwas_stat/'

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

plink.gwas <- subset.plink('1KG_EUR/chr' %&&% chr.input,
                           chr.input, ld.lb.input, ld.ub.input, temp.dir)

if(nrow(eqtl.tab) == 0) {
    write_tsv(data.frame(), path = out.pat.file)
    write_tsv(data.frame(), path = out.mat.file)
    log.msg('Empty QTL file\n')
    q()
}

.run <- function(.gwas.file, .out.file, .gwas) {
    if(!file.exists(.out.file)) {

        .gwas.tab <- read.gwas(.gwas.file)

        .matched <- .gwas.tab %>%
            match.allele(plink.obj = plink.gwas, qtl.tab = eqtl.tab) %>%
                mutate(qtl.beta = qtl.beta / pmax(qtl.se, 1e-8), qtl.se = 1)

        .data <- .matched %>%
            make.zqtl.data(n.permuted = 30)
        gc()

        if(is.null(.data)) {
            write_tsv(data.frame(), path = .out.file)
            log.msg('Empty data\n')
            return(NULL)
        }

        vb.opt <- list(pi = -1, tau = -5,
                       do.hyper = FALSE, tol = 1e-8,
                       gammax = gammax.input, nsingle = 100,
                       vbiter = 5000, do.stdize = TRUE, eigen.tol = eig.tol,
                       rate = 1e-2, nsample = 10, print.interv = 500,
                       weight = FALSE, do.rescale = TRUE,
                       multivar.mediator = TRUE)
        
        z.out <- .data %>%
            run.cammel(xx.gwas = plink.gwas$BED, xx.med = plink.gwas$BED, opt = vb.opt)
        
        var.tab <- get.var.tab(z.out$var.decomp, .data$mediators) %>%
            select(med.id, var.mediated, var.direct.tot)
        
        summary.tab <- .matched %>% group_by(med.id) %>% slice(which.max(abs(qtl.z))) %>%
            select(med.id, gwas.p, gwas.z)
        
        out.tab <- melt.effect(z.out$param.mediated, .data$mediators, .gwas) %>%
            rename(med.id = Var1, gwas = Var2) %>%
                left_join(var.tab) %>%
                    left_join(summary.tab)
        
        out.tab <- out.tab %>%
            mutate(ld.idx = ld.idx,
                   gwas.p.ld = min(.matched$gwas.p),
                   num.genes.ld = nrow(summary.tab))
        
        write_tsv(out.tab, path = .out.file)
    }
}

gwas.file <- gwas.dir %&&% 'ukbb_ad_mat_' %&&% ld.idx %&&% '.txt.gz'
out.file <- out.hdr %&&% '.ad_mat.gz'
.run(gwas.file, out.file, 'ukbb.mat')

gwas.file <- gwas.dir %&&% 'ukbb_ad_pat_' %&&% ld.idx %&&% '.txt.gz'
out.file <- out.hdr %&&% '.ad_pat.gz'
.run(gwas.file, out.file, 'ukbb.pat')

gwas.file <- gwas.dir %&&% 'ukbb_neuroticism_' %&&% ld.idx %&&% '.txt.gz'
out.file <- out.hdr %&&% '.neuroticism.gz'
.run(gwas.file, out.file, 'ukbb.neuroticism')

gwas.file <- gwas.dir %&&% 'ukbb_moodswings_' %&&% ld.idx %&&% '.txt.gz'
out.file <- out.hdr %&&% '.moodswings.gz'
.run(gwas.file, out.file, 'ukbb.moodswings')


log.msg('Successfully finished!\n')
