#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-cammel.R')

if(length(argv) < 4) {
    q()
}

chr.input <- as.integer(argv[1])   # chr.input = 19
ld.lb.input <- as.integer(argv[2]) # ld.lb.input = 8347513
ld.ub.input <- as.integer(argv[3]) # ld.ub.input = 9238393
out.file <- argv[4]                # out.file = 'temp.txt.gz'

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

if(file.exists(out.file)) {
    log.msg('All output files already exist: %s\n', out.file)
    q()
}

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

cis.dist <- 1e6
ctrl.degree <- 3
p.val.cutoff <- 1e-4

expr.dir <- 'ROSMAP_EXPR/count'
geno.dir <- 'ROSMAP_GENO'
covar.file <- 'ROSMAP_EXPR/phenotypes.txt.gz'

library(fqtl)
library(dplyr)
library(tidyr)
library(readr)
library(methods)

gene.info.files <- expr.dir %&&% '/chr' %&&% 1:22 %&&% '.info.gz'
expr.files <- expr.dir %&&% '/chr' %&&% 1:22 %&&% '.data.gz'

.read.tsv <- function(...) suppressMessages(read_tsv(...))

################################################################
## Find potential regulatory genes
info.cols <- c('ensg', 'chr', 'tss', 'tes', 'hgnc')
genes.info <- bind_rows(lapply(gene.info.files, .read.tsv, col_names = info.cols)) %>%
    mutate(hgnc = if_else(is.na(hgnc) | nchar(hgnc) < 1, ensg, hgnc))

################################################################
## Read the whole data and size adjustment using geometric mean
Y <- bind_cols(lapply(expr.files, .read.tsv, col_names = FALSE)) %>%
    adjust.size.factor() %>%
        as.matrix()

covar.tab <- .read.tsv(covar.file)

## remove genes with too many zeros or NAs
.temp <- apply(Y < 1, 2, mean, na.rm = TRUE)
.rm <- which(is.na(.temp) | .temp > .1)

if(length(.rm) > 0) {
    Y <- (Y %c% - .rm) %>% adjust.size.factor()
    genes.info <- (genes.info %r% -.rm) %>%
        mutate(idx = 1:n())
}

################################################################
ld.genes <-
    genes.info %>%
        filter(chr == chr.input) %>%
            filter((tss > (ld.lb.input - cis.dist) & tss < (ld.ub.input + cis.dist)) | (tes > (ld.lb.input - cis.dist) & tes < (ld.ub.input + cis.dist)))

## Include only coding genes
coding.genes <- read_tsv('coding.genes.txt.gz') %>% na.omit()

ld.genes <- ld.genes %>%    
    filter(ensg %in% coding.genes$ensg)

if(nrow(ld.genes) < 1) {
    write_tsv(x = data.frame(), path = out.file)
    log.msg('No gene\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

gene.idx <- ld.genes$idx


################################################################
## Find correlated out-of-chromosome genes to use as control
Y1 <- Y %c% gene.idx %>% stdize.count()
genes.Y1 <- genes.info %r% gene.idx %>%
    mutate(y.col = 1:n()) %>%
        mutate(chr = as.integer(chr)) %>%
            rename(med.id = ensg)

if(nrow(genes.Y1) < 1){
    write_tsv(x = data.frame(), path = out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

################################################################
other.genes <- genes.info %>%
    dplyr::filter(!(chr %in% genes.Y1$chr))

other.idx <- other.genes$idx

Y0 <- Y %c% other.idx

nn.ctrl <- min(ctrl.degree, ncol(Y0))

log.msg('Find %d control genes for each target gene\n', nn.ctrl)

trans.normal <- function(...){
    ret <- log2(1/2 + ...) %>% scale() ## Voom type of transformation
    return(ret)
}

ctrl.idx <- find.cor.idx(Y1 %>% trans.normal(),
                         Y0 %>% trans.normal(),
                         n.ctrl = nn.ctrl,
                         p.val.cutoff = 1e-4)

Y0.ctrl <- Y0 %c% ctrl.idx %>% as.matrix()
genes.ctrl <- other.genes %r% ctrl.idx

################################################################
## Take reference genotypes
plink.hdr <- geno.dir %&&% '/chr' %&&% chr.input

plink.eqtl <- subset.plink(plink.hdr, chr.input, ld.lb.input, ld.ub.input, temp.dir)

if(is.null(plink.eqtl)){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

plink.gwas <- subset.plink('1KG_EUR/chr' %&&% chr.input,
                           chr.input, ld.lb.input, ld.ub.input, temp.dir)

## Read and match two PLINK filesets
plink.matched <- match.plink(plink.gwas, plink.eqtl)

if(is.null(plink.matched)) {
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

plink.gwas <-  plink.matched$gwas
plink.eqtl <-  plink.matched$qtl

if(is.null(plink.eqtl)){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

################################################################
x.fam <- plink.eqtl$FAM %>%
    mutate(x.pos = 1:n()) %>%
        rename(IID = iid)

x.bim <- plink.eqtl$BIM %>%
    mutate(x.col = 1:n())

covar.matched <- covar.tab %>%
    mutate(y.pos = 1:n()) %>%
        left_join(x.fam) %>%
            filter(!is.na(x.pos), !is.na(y.pos))

x.pos <- covar.matched$x.pos
y.pos <- covar.matched$y.pos
pos.df <- data.frame(x.pos, y.pos)

################################################################
xx.std <- plink.eqtl$BED %r% x.pos %>% scale() %>% rm.na.zero()

Y0.std <- Y0.ctrl %>% trans.normal() %r% y.pos %>% scale() %>% rm.na.zero()
colnames(Y0.std) <- 1:ncol(Y0.ctrl)

x.y0.qtl <- calc.qtl.stat(xx.std, Y0.std)

max.qtl <- x.y0.qtl %>%
    rename(gene = y.col) %>%
        group_by(gene)%>%
            slice(which.min(p.val))

rm.cols <- max.qtl %>%
    dplyr::filter(p.val < p.val.cutoff) %>%
        dplyr::select(gene) %>% unique()

if(nrow(rm.cols) > 0){
    rm.idx <- rm.cols$gene
    genes.ctrl <- genes.ctrl %r% (-rm.idx)
    Y0.ctrl <- Y0.ctrl %c% (-rm.idx)
}

log.msg('Identified %d Y0 covariates; removed %d genes\n',
        ncol(Y0.ctrl),
        nrow(rm.cols))

################################################################
covar.mat <- covar.tab %>%
    select(amyloid_sqrt, tangles_sqrt, cog_ep_random_slope, msex, studyn, age_death, educ, RIN.by.BSP) %>%
        as.matrix() %>%
            scale() %>%
                as.data.frame() %>%
                    mutate(intercept = 1) %>%
                        as.matrix()

################################################################
K <- min(10, min(ncol(Y0.ctrl), ncol(covar.mat)))
lodds.cutoff <- log(0.1) - log(0.9)

opt.reg <- list(vbiter = 3500, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi.ub = 0, pi.lb = lodds.cutoff, tau = -4, do.hyper = TRUE, jitter = 1e-2,
                model = 'nb', out.residual = TRUE, k = K,
                svd.init = TRUE, print.interv = 100)

################################################################
## Control genes are available
if(ncol(Y0.ctrl) > 0) {

    ## 0. remove covariance effects (and convert NB to Gaussian)
    y0.out <- fqtl.regress(y = Y0.ctrl,
                           x.mean = covar.mat,
                           factored = TRUE,
                           options = opt.reg)

    y0.covar <- y0.out$resid$theta %>% scale()

} else {
    y0.covar <- matrix(nrow = nrow(covar.mat), ncol = 0)
}

################################################################
## Estimate multivariate models
covar.mat.combined <- cbind(covar.mat, y0.covar)

xx.std <- plink.eqtl$BED %r% pos.df$x.pos %>% scale()

opt.reg <- list(vbiter = 3500, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi.ub = 0, pi.lb = lodds.cutoff, tau = -4, do.hyper = TRUE, jitter = 1e-2,
                model = 'nb', out.residual = FALSE, print.interv = 100)

y1.out <- fqtl.regress(y = Y1 %r% pos.df$y.pos,
                       x.mean = xx.std,
                       c.mean = covar.mat.combined %r% pos.df$y.pos,
                       opt = opt.reg)

out.tab <- effect2tab(y1.out$mean) %>%
    left_join(x.bim) %>%
        left_join(genes.Y1) %>%
            rename(qtl.beta = theta, qtl.se = theta.se, qtl.lodds = lodds) %>%
                rename(qtl.a1 = plink.a1, qtl.a2 = plink.a2) %>%
                    select(chr, rs, snp.loc, med.id, dplyr::starts_with('qtl'))

valid.med <- out.tab %>%
    group_by(med.id) %>%
        slice(which.max(qtl.lodds)) %>%
            filter(qtl.lodds > lodds.cutoff) %>%
                select(med.id) %>%
                    .unlist()

qtl.out <- out.tab %>%
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
