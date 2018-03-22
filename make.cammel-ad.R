#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

if(length(argv) < 7) {
    q()
}

ld.idx <- as.integer(argv[1])           # e.g., ld.idx = 133
qtl.file <- argv[2]                     # e.g., qtl.file = 'cis-eqtl/geuvadis/133_qtl.txt.gz'
qtl.sample.size <- as.integer(argv[3])  # e.g., qtl.sample.size = 500
geno.dir <- argv[4]                     # e.g., geno.dir = '1KG_EUR' # (eQTL genotype matrix)
gammax.input <- as.numeric(argv[5])     # e.g., gammax.input = 1e4, 1e3, 1e2
eig.tol <- as.numeric(argv[6])          # e.g., eig.tol = 1e-1
out.hdr <- argv[7]                      # e.g., out.hdr = 'temp.cammel'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

library(zqtl)
library(dplyr)
library(readr)
library(methods)

################################################################
igap.sample.size <- 74056
n.null <- 3

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

################################################################

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

## Read QTL statistics and measure basic statistics
## chr rs snp.loc med.id qtl.a1 qtl.a2 qtl.beta qtl.z
## i   c  i       c      c      c      d        d
qtl.tab <- read_tsv(qtl.file, col_types = 'icicccdd')

if(nrow(qtl.tab) == 0) {
    write_tsv(data.frame(), path = out.tab.file)
    write_tsv(data.frame(), path = out.null.file)
    log.msg('Empty QTL file\n')
    q()
}

qtl.tab <- qtl.tab %>% dplyr::select(-rs)

med.qtl.stat <- qtl.tab %>%
    group_by(med.id) %>%
        summarize(n.qtl.4 = sum(abs(qtl.z) > 4))

igap.gwas.tab <- read_tsv(igap.gwas.file)

################################################################
igap.matched <- igap.gwas.tab %>%
    match.allele(plink.obj = plink.eqtl, qtl.tab = qtl.tab)

igap.data <- igap.matched %>%
    make.zqtl.data()

################################################################
run.cammel.null <- function(xx.gwas, xx.med, zqtl.data, n.null, gwas.sample.size) {

    vb.opt <- list(pi.ub = -1, pi.lb = -5, tau = -5, do.hyper = TRUE, tol = 1e-8, gammax = gammax.input,
                   vbiter = 3000, do.stdize = TRUE, eigen.tol = eig.tol,
                   rate = 1e-2, nsample = 10, print.interv = 1500,
                   weight = FALSE, do.rescale = FALSE)

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()
    ret <- NULL
    for(r in 1:n.null) {

        gwas.null <- make.zqtl.null(xx.gwas.mat, zqtl.data$gwas.se, eig.tol = eig.tol)
        qtl.null <- make.zqtl.null(xx.med.mat, zqtl.data$qtl.se, eig.tol = eig.tol)

        z.out <- fit.med.zqtl(gwas.null, zqtl.data$gwas.se,
                              qtl.null, zqtl.data$qtl.se,
                              X.gwas = xx.gwas.mat, X.med = xx.med.mat,
                              n = gwas.sample.size, n.med = qtl.sample.size,
                              options = vb.opt)

        null.stat <- melt.effect(z.out$param.mediated, zqtl.data$mediators, r) %>%
            mutate(theta.var = sqrt(theta.var)) %>% rename(theta.se = theta.var) %>%
                mutate(log.var = as.numeric(signif(log(z.out$var.decomp$var.med.each), 2)))

        print(null.stat %>% filter(lodds > 0) %>% as.data.frame())
        ret <- bind_rows(ret, null.stat)
        log.msg('\nnull round = %d / %d\n', r, n.null)
    }
    ret <- ret %>% rename(med.id = Var1, null = Var2)
    return(ret)
}

run.cammel <- function(xx.gwas, xx.med, zqtl.data, gwas.sample.size) {

    vb.opt <- list(pi.ub = -1, pi.lb = -5, tau = -5, do.hyper = TRUE, tol = 1e-8, gammax = gammax.input,
                   vbiter = 3000, do.stdize = TRUE, eigen.tol = eig.tol,
                   rate = 1e-2, nsample = 10, print.interv = 1500,
                   weight = FALSE, do.rescale = FALSE)

    if(is.null(zqtl.data)) return(NULL)

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()

    z.out <- fit.med.zqtl(zqtl.data$gwas.beta, zqtl.data$gwas.se,
                          zqtl.data$qtl.beta, zqtl.data$qtl.se,
                          X.gwas = xx.gwas.mat,
                          X.med = xx.med.mat,
                          n = gwas.sample.size,
                          n.med = qtl.sample.size, options = vb.opt)
    return(z.out)
}

get.var.tab <- function(var.decomp, mediators) {

    if(is.null(var.decomp)) return(NULL)

    ret <- data.frame(med.id = mediators,
                      var.mediated = signif(var.decomp$var.med.each, 2),
                      var.mediated.tot = signif(var.decomp$var.med.mean, 2),
                      var.mediated.tot.se = signif(sqrt(var.decomp$var.med.var), 2),
                      var.direct.tot = signif(var.decomp$var.direct.mean, 2),
                      var.direct.tot.se = signif(sqrt(var.decomp$var.direct.var), 2))
    return(ret)
}

## gene-level QTL and GWAS stat summary
get.summary.tab <- function(.gwas.tab) {
    if(is.null(.gwas.tab)) return(NULL)
    .temp <- .gwas.tab %>%
        match.allele(plink.obj = plink.eqtl, qtl.tab = qtl.tab)

    .temp.gwas <- .temp %>% group_by(med.id) %>%
        slice(which.min(gwas.p)) %>%
            select(med.id, rs, snp.loc, gwas.p, gwas.beta, gwas.se, qtl.z)

    .temp.qtl <- .temp %>% group_by(med.id) %>%
        slice(which.max(qtl.z)) %>%
            select(med.id, rs, snp.loc, gwas.p, gwas.beta, gwas.se, qtl.z)

    ret <- .temp.gwas %>%
        left_join(.temp.qtl, by = 'med.id', suffix = c('.by.gwas', '.by.qtl'))

    ret <- ret %>% dplyr::select_(.dots = sort(names(ret)))

    return(ret)
}

get.effect.tab <- function(z.out, z.data, gwas.tab, data.name) {
    if(is.null(z.out)) return(NULL)
    z.effect <- melt.effect(z.out$param.mediated, z.data$mediators, data.name) %>%
        rename(med.id = Var1, gwas = Var2) %>% left_join(med.qtl.stat) %>%
            left_join(get.var.tab(z.out$var.decomp, z.data$mediators)) %>%
                left_join(get.summary.tab(gwas.tab))
    return(z.effect)
}

gc()

igap.out <- run.cammel(plink.gwas$BED, plink.eqtl$BED, igap.data, igap.sample.size)
igap.null <- run.cammel.null(plink.gwas$BED, plink.eqtl$BED, igap.data, n.null, igap.sample.size)

################################################################
## collect effect sizes
igap.effect <- get.effect.tab(igap.out, igap.data, igap.gwas.tab, 'igap')

out.tab <- igap.effect %>%
    dplyr::mutate(chr = chr.input, ld.lb = ld.lb.input, ld.ub = ld.ub.input,
                  n.snp = nrow(plink.eqtl$BIM))

null.tab <- igap.null %>% mutate(gwas = 'igap')

write_tsv(out.tab, path = out.tab.file)
write_tsv(null.tab, path = out.null.file)

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
