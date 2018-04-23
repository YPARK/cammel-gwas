#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')

if(length(argv) != 7) {
    q()
}

src.chr.input <- as.integer(argv[1])   # src.chr.input = 19
src.ld.lb.input <- as.integer(argv[2]) # src.ld.lb.input = 8347513
src.ld.ub.input <- as.integer(argv[3]) # src.ld.ub.input = 9238393

tgt.chr.input <- as.integer(argv[4])   # tgt.chr.input = 1
tgt.ld.lb.input <- as.integer(argv[5]) # tgt.ld.lb.input = 10583
tgt.ld.ub.input <- as.integer(argv[6]) # tgt.ld.ub.input = 1892607

out.file <- argv[7] # out.file = 'temp.txt.gz'

######################
## source -> target ##
######################

################################################################
temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel.trans/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel.trans/' %&&% out.file %&&%
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
library(readr)
library(methods)

gene.info.files <- expr.dir %&&% '/chr' %&&% 1:22 %&&% '.info.gz'
expr.files <- expr.dir %&&% '/chr' %&&% 1:22 %&&% '.data.gz'

################################################################
## Find potential regulatory genes
info.cols <- c('ensg', 'chr', 'tss', 'tes', 'hgnc')
genes.info <- bind_rows(lapply(gene.info.files, read_tsv, col_names = info.cols)) %>%
    mutate(hgnc = if_else(is.na(hgnc) | nchar(hgnc) < 1, ensg, hgnc))


################################################################
## Read the whole data and size adjustment using geometric mean
Y <- lapply(expr.files, read_tsv, col_names = FALSE) %>%
    bind_cols() %>%
        adjust.size.factor() %>%
            as.matrix()

covar.tab <- read_tsv(covar.file)

## remove genes with too many zeros or NAs
.temp <- apply(Y < 1, 2, mean, na.rm = TRUE)
.rm <- which(is.na(.temp) | .temp > .1)

if(length(.rm) > 0) {
    Y <- (Y %c% - .rm) %>% adjust.size.factor()
    
    genes.info <- (genes.info %r% -.rm) %>%
        mutate(idx = 1:n())
}

################################################################
## 1. take genes within source and target regions

take.genes <- function(.chr, .lb, .ub, .dist = cis.dist) {
    genes.info %>% dplyr::filter(chr == .chr) %>%
        dplyr::filter(tss > .lb - .dist, tes < .ub + .dist)
}

src.genes <- take.genes(src.chr.input, src.ld.lb.input, src.ld.ub.input)
tgt.genes <- take.genes(tgt.chr.input, tgt.ld.lb.input, tgt.ld.ub.input)


if(nrow(tgt.genes) < 1) {
    write_tsv(x = data.frame(), path = out.file)
    log.msg('No gene\n')
    system('rm -r ' %&&% temp.dir)
    q()
}

################################################################
## 2. find correlated genes with gene expressions of the source and
## target regions

Y.src <- Y %c% src.genes$idx
Y.tgt <- Y %c% tgt.genes$idx

chr.taboo <- unique(c(src.genes$chr, tgt.genes$chr))

## something out of the chromosomes
other.genes <- genes.info %>% dplyr::filter(!(chr %in% chr.taboo))

Y0 <- Y %c% other.genes$idx

nn.ctrl <- min(ctrl.degree, ncol(Y0))

trans.normal <- function(...){
    ret <- log2(1/2 + ...) %>% scale() ## Voom type of transformation
    return(ret)
}

ctrl.src.idx <- find.cor.idx(Y.src %>% trans.normal(),
                             Y0 %>% trans.normal(),
                             n.ctrl = nn.ctrl,
                             p.val.cutoff = 1e-4)

ctrl.tgt.idx <- find.cor.idx(Y.tgt %>% trans.normal(),
                             Y0 %>% trans.normal(),
                             n.ctrl = nn.ctrl,
                             p.val.cutoff = 1e-4)

## they are on the same Y0
ctrl.idx <- unique(c(ctrl.src.idx, ctrl.tgt.idx))

Y0.ctrl <- Y0 %c% ctrl.idx %>% as.matrix()
genes.ctrl <- other.genes %r% ctrl.idx

################################################################
## 3. But filter out genes correlated with genotype information within
## the source region.

plink.hdr <- geno.dir %&&% '/chr' %&&% src.chr.input

plink <- subset.plink(plink.hdr, src.chr.input, src.ld.lb.input, src.ld.ub.input, temp.dir)

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

covar.matched <- covar.tab %>%
    mutate(y.pos = 1:n()) %>%
        left_join(x.fam) %>%
            filter(!is.na(x.pos), !is.na(y.pos))

x.pos <- covar.matched$x.pos
y.pos <- covar.matched$y.pos
pos.df <- data.frame(x.pos, y.pos)

################################################################
xx.std <- plink$BED %r% x.pos %>% scale() %>% rm.na.zero()

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
opt.reg <- list(vbiter = 5000, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi = -1, tau = -4, do.hyper = FALSE, jitter = 0.01,
                model = 'nb', out.residual = TRUE, k = 10,
                svd.init = TRUE, print.interv = 100)

## FQTL if there is suffcient sharing exists
Y1 <- Y.tgt %>% stdize.count()

genes.Y1 <- tgt.genes %>%
    mutate(y.col = 1:n()) %>%
        mutate(chr = as.integer(chr)) %>%
            rename(med.id = ensg)

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
            left_join(genes.Y1, by = 'y.col', suffix = c('.snp', '.gene')) %>%
                dplyr::rename(qtl.a1 = plink.a1, qtl.a2 = plink.a2) %>%
                    mutate(qtl.beta = signif(beta, 4), qtl.z = signif(beta/se, 4))

    ret <- ret %>%
        dplyr::select(chr.snp, rs, snp.loc, chr.gene, med.id, dplyr::starts_with('qtl'))

    ret <- ret %>%
        dplyr::filter(nchar(qtl.a1) == 1, nchar(qtl.a2) == 1)

    return(ret)
}

################################################################
## Control genes are available
if(ncol(Y0.ctrl) > 0) {

    ## remove covariance effects (and convert NB to Gaussian)
    y0.out <- fqtl.regress(y = Y0.ctrl,
                           x.mean = covar.mat,
                           options = opt.reg)

    y0.covar <- y0.out$resid$theta %>% scale()

    ## report QTL statistics
    y1.out <- fqtl.regress(y = Y1 %r% pos.df$y.pos,
                           x.mean = y0.covar %r% pos.df$y.pos,
                           c.mean = covar.mat %r% pos.df$y.pos,
                           factored = .factored,
                           options = opt.reg)

    yy.std <- y1.out$resid$theta %>% scale()
    xx.std <- plink$BED %r% pos.df$x.pos %>% scale()
    overall.qtl <- .get.marginal.qtl(xx.std, yy.std)

} else {

    ## overall samples
    y1.out <- fqtl.regress(y = Y1 %r% pos.df$y.pos,
                           x.mean = covar.mat %r% pos.df$y.pos,
                           factored = .factored,
                           options = opt.reg)

    ## report QTL statistics
    yy.std <- y1.out$resid$theta %>% scale()
    xx.std <- plink$BED %r% pos.df$x.pos %>% scale()
    overall.qtl <- .get.marginal.qtl(xx.std, yy.std)

}

write_tsv(x = overall.qtl, path = out.file)

system('rm -r ' %&&% temp.dir)
log.msg('Successfully finished!\n')
