#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-cammel.R')

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
            filter((tss > (ld.lb.input - cis.dist) & tss < (ld.ub.input + cis.dist)) | (tes > (ld.lb.input - cis.dist) & tes < (ld.ub.input + cis.dist)))

## Include only coding genes
coding.genes <- read_tsv('coding.genes.txt.gz') %>% na.omit()

ld.genes <- ld.genes %>%    
    filter(med.id %in% coding.genes$ensg)

if(nrow(ld.genes) < 1) {
    log.msg('No valid gene\n')
    write_tsv(data.frame(), path = out.file)
    q()
}

gene.idx <- ld.genes$idx

################################################################
## take genotypes

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel.mayo/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel.mayo/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

plink.eqtl <- subset.plink(plink.hdr, chr.input, ld.lb.input, ld.ub.input, temp.dir)

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
                         p.val.cutoff = 1e-4)

Y0.ctrl <- Y0 %c% ctrl.idx %>% as.matrix() %>% stdize.count()
genes.ctrl <- other.genes %r% ctrl.idx

x.pos <- samples.info$plink.idx

################################################################
## check genetic correlation
xx.std <- plink.eqtl$BED %r% x.pos %>% scale() %>% rm.na.zero()
colnames(xx.std) <- plink.eqtl$BIM$snp.loc
x.bim <- plink.eqtl$BIM %>%
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

xx.std <- plink.eqtl$BED %r% x.pos %>% scale()

opt.reg <- list(vbiter = 3500, gammax = 1e4, tol = 1e-8, rate = 1e-2,
                pi.ub = 0, pi.lb = lodds.cutoff, tau = -4, do.hyper = TRUE, jitter = 1e-2,
                model = 'nb', out.residual = FALSE, print.interv = 100)

y1.out <- fqtl.regress(y = Y1,
                       x.mean = xx.std,
                       c.mean = covar.mat.combined,
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
