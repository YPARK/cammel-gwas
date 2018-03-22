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

expr.file <- 'cis-eqtl/geuvadis_data/expression.txt.gz'
sample.file <- 'cis-eqtl/geuvadis_data/samples.txt.gz'
gene.file <- 'cis-eqtl/geuvadis_data/genes.txt.gz'
geno.dir <- '1KG_EUR'
covar.file <- '1KG_EUR/1kg_sample_info.txt.gz'

cis.dist <- 1e6
ctrl.degree <- 3
p.val.cutoff <- 1e-4

library(fqtl)
library(dplyr)
library(readr)
library(methods)

covar.tab.full <- read_tsv(covar.file)
genes.info <- read_tsv(gene.file) %>% mutate(idx = 1:n())
sample.info <- read_tsv(sample.file)

Y <- read_tsv(expr.file) %>% as.matrix() %>% t()

################################################################
ld.genes <-
    genes.info %>%
        filter(chr == chr.input) %>%
            filter(gene.loc > (ld.lb.input - cis.dist), gene.loc < (ld.ub.input + cis.dist))

gene.idx <- ld.genes$idx

if(nrow(ld.genes) < 1) {
    log.msg('No valid gene\n')
    write_tsv(data.frame(), path = out.file)
    q()
}

################################################################
## Only consider European samples
covar.tab <-
    sample.info %>%
        left_join(covar.tab.full) %>%
            dplyr::rename(IID = Sample)

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel.geuvadis/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel.geuvadis/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

plink.hdr <- geno.dir %&&% '/chr' %&&% chr.input

plink <- subset.plink(plink.hdr, chr.input, ld.lb.input, ld.ub.input, temp.dir)

if(is.null(plink)){
    write_tsv(data.frame(), path = out.file)
    log.msg('Just wrote empty QTL files!\n')
    system('rm -r ' %&&% temp.dir)    
    q()
}

x.fam <- plink$FAM %>%
    mutate(x.pos = 1:n()) %>%
        rename(IID = iid)

x.bim <- plink$BIM %>%
    mutate(x.col = 1:n())

.matched <- covar.tab %>%
    mutate(y.pos = 1:n()) %>%
        left_join(x.fam) %>%
            filter(!is.na(x.pos), !is.na(y.pos))

.temp <- apply(.matched[, -1], 2, sd)
rm.cols <- which(is.na(.temp) | .temp < 1e-4)

x.pos <- .matched$x.pos
y.pos <- .matched$y.pos
pos.df <- data.frame(x.pos, y.pos)

covar.matched <- .matched %c% -(rm.cols + 1)

################################################################
## Find correlated out-of-chromosome genes to use as control
Y1 <- Y %c% gene.idx %r% pos.df$y.pos %>% stdize.count()

genes.Y1 <- genes.info %r% gene.idx %>% mutate(y.col = 1:n()) %>%
    mutate(chr = as.integer(chr)) %>%
        rename(med.id = ensg)

if(nrow(genes.Y1) < 1){
    write_tsv(x = data.frame(), path = out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

################################################################
other.genes <- genes.info %>%
    dplyr::filter(chr != chr.input, chr != 'X', chr != 'Y')

other.idx <- other.genes$idx

Y0 <- Y %c% other.idx %r% pos.df$y.pos

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

################################################################
xx.std <- plink$BED %r% x.pos %>% scale() %>% rm.na.zero()

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
    rm.idx <- rm.cols$gene
    genes.ctrl <- genes.ctrl %r% (-rm.idx)
    Y0.ctrl <- Y0.ctrl %c% (-rm.idx)
}

log.msg('Identified %d Y0 covariates; removed %d genes\n',
        ncol(Y0.ctrl),
        nrow(rm.cols))

covar.mat <- covar.matched %>%
    select(-x.pos, -y.pos, -IID, -n) %>%
        scale() %>%
            as.data.frame() %>%
                mutate(intercept = 1) %>%
                    as.matrix()

################################################################
opt.reg <- list(vbiter = 5000, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi = -1, tau = -4, do.hyper = FALSE, jitter = 0.01,
                model = 'nb', out.residual = TRUE, k = 10,
                svd.init = TRUE, print.interv = 100)

## FQTL if there is suffcient sharing exists
if(ncol(Y1) > opt.reg$k) {
    .factored <- TRUE
    if(ncol(Y0.ctrl) > opt.reg$k) {
        opt.reg$svd.init <- TRUE
        log.msg('Will use SVD initialization\n')
    } else {
        opt.reg$svd.init <- FALSE
        log.msg('Will use random initialization\n')
    }
} else {
    .factored <- FALSE
    log.msg('Will not fit FQTL\n')
}

.get.marginal.qtl <- function(...) {
    ret <- calc.qtl.stat(...) %>%
        left_join(x.bim) %>%
            left_join(genes.Y1) %>%
                rename(qtl.a1 = plink.a1, qtl.a2 = plink.a2) %>%
                    mutate(qtl.beta = signif(beta, 4), qtl.z = signif(beta/se, 4)) %>%
                        select(chr, rs, snp.loc, med.id, dplyr::starts_with('qtl'))
    return(ret)
}

################################################################
## Control genes are available
if(ncol(Y0.ctrl) > 0) {

    ## 0. remove covariance effects (and convert NB to Gaussian)
    y0.out <- fqtl.regress(y = Y0.ctrl,
                           x.mean = covar.mat,
                           options = opt.reg)

    y0.covar <- y0.out$resid$theta %>% scale()

    ## Correction
    y1.out <- fqtl.regress(y = Y1, 
                           x.mean = y0.covar,
                           c.mean = covar.mat,
                           factored = .factored,
                           options = opt.reg)

    ## Report QTL statistics
    yy.std <- y1.out$resid$theta %>% scale()
    xx.std <- plink$BED %r% pos.df$x.pos %>% scale()
    qtl.out <- .get.marginal.qtl(xx.std, yy.std)

} else {

    ## Control genes are NOT available
    ## Correction
    y1.out <- fqtl.regress(y = Y1, 
                           x.mean = covar.mat,
                           factored = .factored,
                           options = opt.reg)

    ## Report QTL statistics
    yy.std <- y1.out$resid$theta %>% scale()
    xx.std <- plink$BED %r% pos.df$x.pos %>% scale()
    qtl.out <- .get.marginal.qtl(xx.std, yy.std)

}

write_tsv(x = qtl.out, path = out.file)

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
