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
out.file <- argv[4]                # out.file = 'temp.txt.gz'

if(file.exists(out.file)) {
    log.msg('File exists: %s\n', out.file)
    q()
}

cis.dist <- 1e6
ctrl.degree <- 3
p.val.cutoff <- 1e-4

library(fqtl)
library(dplyr)
library(tidyr)
library(readr)
library(methods)

expr.file <- 'cis-eqtl/mayo_data/expression.txt.gz'
sample.file <- 'cis-eqtl/mayo_data/samples.txt.gz'
gene.file <- 'cis-eqtl/mayo_data/genes.txt.gz'
plink.hdr <- 'Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed'

expr.data <- read_tsv(expr.file) %>% t() %>% as.matrix()
samples.info <- read_tsv(sample.file)
genes.info <- read_tsv(gene.file) %>%
    mutate(idx = 1:n()) %>%
        rename(med.id = ensembl_id)

covar.tab <- samples.info %>%
    select(ID, RIN, Diagnosis, Gender, AgeAtDeath, ApoE, PMI) %>%
        mutate(Gender = if_else(Gender == 'M', 1, -1)) %>%
            mutate(AgeAtDeath = if_else(AgeAtDeath == '90_or_above', '90', AgeAtDeath)) %>%
                mutate(AgeAtDeath = as.numeric(AgeAtDeath)) %>%
                    mutate(val = 1) %>%
                        spread(key = Diagnosis, value = val, fill = 0) %>%
                            mutate(val = 1) %>%
                                spread(key = ApoE, value = val, fill = 0)

.order <- match(samples.info$ID, covar.tab$ID)

covar.mat <- covar.tab %r% .order %>% select(-ID) %>%
    scale() %>%
        as.data.frame() %>%
            mutate(intercept = 1) %>%
                as.matrix()

################################################################
ld.genes <-
    genes.info %>%
        filter(chr == chr.input) %>%
            filter(tss > (ld.lb.input - cis.dist), tes < (ld.ub.input + cis.dist))

gene.idx <- ld.genes$idx

if(nrow(ld.genes) < 1) {
    log.msg('No valid gene\n')
    write_tsv(data.frame(), path = out.file)
    q()
}

################################################################
## take genotypes

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel.mayo/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel.mayo/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

plink <- subset.plink(plink.hdr, chr.input, ld.lb.input, ld.ub.input, temp.dir)

if(is.null(plink)){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)    
    q()
}

################################################################
## Find correlated out-of-chromosome genes to use as control
Y1 <- expr.data %c% gene.idx %>% stdize.count()
genes.Y1 <- genes.info %r% gene.idx %>% mutate(y.col = 1:n())

if(nrow(genes.Y1) < 1){
    write_tsv(x = data.frame(), path = out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

other.genes <- genes.info %>%
    dplyr::filter(chr != chr.input, chr != 'X', chr != 'Y')

other.idx <- other.genes$idx

Y0 <- expr.data %c% other.idx

nn.ctrl <- min(ctrl.degree, ncol(Y0))

log.msg('Find %d control genes for each target gene\n', nn.ctrl)

trans.normal <- function(...){
    ret <- log2(1/2 + ...) %>% scale() ## Voom type of transformation
    return(ret)
}

ctrl.idx <- find.cor.idx(Y1 %>% trans.normal(),
                         Y0 %>% trans.normal(),
                         n.ctrl = nn.ctrl,
                         p.val.cutoff = 1e-2)

Y0.ctrl <- Y0 %c% ctrl.idx %>% as.matrix() %>% stdize.count()
genes.ctrl <- other.genes %r% ctrl.idx

x.pos <- samples.info$plink.idx

################################################################
## check genetic correlation
xx.std <- plink$BED %r% x.pos %>% scale() %>% rm.na.zero()
colnames(xx.std) <- plink$BIM$snp.loc
x.bim <- plink$BIM %>%
    mutate(x.col = 1:n())

Y0.std <- Y0.ctrl %>% trans.normal() %>% scale() %>% rm.na.zero()
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
    rm.idx <- rm.cols$gene %>% as.integer()
    genes.ctrl <- genes.ctrl %r% (-rm.idx)
    Y0.ctrl <- Y0.ctrl %c% (-rm.idx)
}

log.msg('Identified %d Y0 covariates; removed %d genes\n',
        ncol(Y0.ctrl),
        nrow(rm.cols))

################################################################
opt.reg <- list(vbiter = 5000, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi = -1, tau = -4, do.hyper = FALSE, jitter = 0.01,
                model = 'nb', out.residual = TRUE, k = 10,
                svd.init = TRUE, print.interv = 100)

################################################################
## Control genes are available
if(ncol(Y0.ctrl) > 0) {

    ## 0. remove covariance effects (and convert NB to Gaussian)
    y0.out <- fqtl.regress(y = Y0.ctrl,
                           x.mean = covar.mat,
                           options = opt.reg)

    y0.covar <- y0.out$resid$theta %>% scale()

} else {
    y0.covar <- matrix(nrow = nrow(covar.mat), ncol = 0)
}

## Estimate multivariate models
covar.mat.combined <- cbind(covar.mat, y0.covar)

xx.std <- plink$BED %r% x.pos %>% scale()

opt.reg <- list(vbiter = 5000, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi = -1, tau = -4, do.hyper = FALSE, jitter = 0.01,
                model = 'nb', out.residual = FALSE, k = 10,
                svd.init = TRUE, print.interv = 100)

y1.out <- fqtl.regress(y = Y1,
                       x.mean = xx.std,
                       c.mean = covar.mat.combined,
                       opt = opt.reg)

resid <- y1.out$resid$theta %>% scale()

pip.cutoff <- 0.9
logit <- function(x) log(x) - log(1 - x)
lodds.cutoff <- logit(pip.cutoff)

mat2tab <- function(mat, val.name) {

    ret <- as.matrix(mat) %>% as.data.frame()
    col.names <- 1:ncol(ret)
    colnames(ret) <- col.names

    ret <- ret %>%
        mutate(x.col = 1:n()) %>%
            gather_(key_col= 'y.col', value_col = val.name, col.names)

    return(ret)
}

effect2tab <- function(effect, lodds.cutoff) {
    ret <- effect$lodds %>% mat2tab(val.name = 'lodds') %>%
        filter(lodds > lodds.cutoff) %>%
            left_join(effect$theta %>% mat2tab(val.name = 'theta')) %>%
                left_join(effect$theta.var %>% mat2tab(val.name = 'theta.var')) %>%
                    mutate(theta.se = sqrt(theta.var)) %>%
                        select(-theta.var) %>%
                            filter(abs(theta.se) > 1e-8) %>%
                                mutate(x.col = as.integer(x.col)) %>%
                                    mutate(y.col = as.integer(y.col))
}

mult.tab <- effect2tab(y1.out$mean, lodds.cutoff) %>%
    left_join(genes.Y1)

if(nrow(mult.tab) < 1){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

X <- plink$BED

multi.to.uni <- function(tab) {
    eta <- (X %c% tab$x.col) %*% matrix(as.numeric(tab$theta), ncol = 1)
    ret <- calc.qtl.stat(X, eta) %>%
        left_join(x.bim) %>%
            rename(qtl.a1 = plink.a1, qtl.a2 = plink.a2) %>%
                mutate(qtl.beta = signif(beta, 4), qtl.z = signif(beta/se, 4))
}

qtl.out <- mult.tab %>% group_by(med.id) %>%
    do(stat = multi.to.uni(.)) %>%
        unnest() %>%
            select(chr, rs, snp.loc, med.id, dplyr::starts_with('qtl'))

write_tsv(x = qtl.out, path = out.file)

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
