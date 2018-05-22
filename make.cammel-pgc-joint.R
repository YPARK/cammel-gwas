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

################################################################
## Run CaMMEL on PGC GWAS statistics
dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

gwas.data <- c('pgc_scz', 'pgc_ocd', 'pgc_adhd', 'pgc_asd', 'pgc_ptsd_civ_ea', 'pgc_mdd')
gwas.dir <- './gwas_stat/'
gwas.files <- gwas.dir %&&% gwas.data %&&% '_' %&&% ld.idx %&&% '.txt.gz'

if(!all(sapply(gwas.files, file.exists))) {
    log.msg('Missing input files: %s\n', paste(gwas.files, collapse = ', '))
    q()
}

out.files <- out.hdr %&&% '.' %&&% gwas.data %&&% '.gz'

if(all(sapply(out.files, file.exists))) {
    log.msg('all the output files exist: %s\n', paste(out.files, collapse = ', '))
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
eqtl.tab <- read.multivar.eqtl(eqtl.data, eqtl.data.files)
log.msg('Read %d rows eQTL tab\n', nrow(eqtl.tab))

################################################################
ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)
chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[ld.idx, 'chr']) %>%
    as.integer()
ld.lb.input <- ld.info.tab[ld.idx, 'start'] %>% as.integer()
ld.ub.input <- ld.info.tab[ld.idx, 'stop'] %>% as.integer()

################################################################
## Read plink
temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel-ad/' %&&% out.hdr %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

plink.gwas <- subset.plink('1KG_EUR/chr' %&&% chr.input,
                           chr.input, ld.lb.input, ld.ub.input, temp.dir)

system('[ -d ' %&&% temp.dir %&&% ' ] && rm -r ' %&&% temp.dir)

.run <- function(.gwas.file, .out.file, .gwas) {
    if(!file.exists(.out.file)) {

        .gwas.tab <- read.gwas(.gwas.file)

        .matched <- .gwas.tab %>%
            match.allele(plink.obj = plink.gwas, qtl.tab = eqtl.tab) %>%
                mutate(qtl.beta = qtl.beta / pmax(qtl.se, 1e-8), qtl.se = 1)
        
        .data <- .matched %>%
            make.zqtl.data(n.permuted = 20)
        gc()

        if(is.null(.data)) {
            write_tsv(data.frame(), path = .out.file)
            log.msg('Empty data\n')
            return(NULL)
        }

        vb.opt <- list(pi.ub = -1/2, pi.lb = -3, tau = -5,
                       do.hyper = TRUE, tol = 1e-8,
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

for(j in 1:length(gwas.files)) {
    .run(gwas.files[j], out.files[j], gwas.data[j])
}
log.msg('Successfully finished!\n')
