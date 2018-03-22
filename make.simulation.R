#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

################################################################
options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-sim.R')

if(length(argv) < 10) { q() }

ld.idx <- as.integer(argv[1])                 ## ld.idx = 133        # LD index
sample.ratio <- as.numeric(.eval(argv[2]))    ## sample.ratio = 0.3  # Ratio between eQTL and GWAS sample size
n.genes <- as.numeric(.eval(argv[3]))         ## n.genes = 100       # Ratio between SNPs and genes
n.causal.genes <- as.integer(.eval(argv[4]))  ## n.causal.genes = 3  # causal genes
n.eqtl.per.gene <- as.integer(.eval(argv[5])) ## n.eqtl.per.gene = 5 # causal eQTLs per gene
n.direct.snp <- as.integer(.eval(argv[6]))    ## n.direct.snp = 1
h2.eqtl <- as.numeric(.eval(argv[7]))         ## h2.eqtl = 0.17      # eQTL heritability
h2.med <- as.numeric(.eval(argv[8]))          ## h2.med = 0.2        # mediation heritability
h2.unmed <- as.numeric(.eval(argv[9]))        ## h2.unmed = 0.3      # unmediated heritability
out.file <- argv[10]                           ## out.file = 'temp.txt.gz'

if(file.exists(out.file)) {
    log.msg('file exists: %s', out.file)
    q()
}

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

library(zqtl)
library(dplyr)
library(tidyr)
library(readr)
library(glmnet)
library(methods)

ld.info.file <- 'ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)

chr.input <- ld.info.tab[ld.idx, 'chr'] %>%
    gsub(pattern = 'chr', replacement = '') %>%
        as.integer()

ld.lb.input <- ld.info.tab[ld.idx, 'start'] %>% as.integer()
ld.ub.input <- ld.info.tab[ld.idx, 'stop'] %>% as.integer()

temp.dir <- system('mkdir -p /broad/hptmp/ypp/cammel-simulation/' %&&% out.file %&&%
                   '; mktemp -d /broad/hptmp/ypp/cammel-simulation/' %&&% out.file %&&%
                   '/temp.XXXXXXXX',
                   intern = TRUE,
                   ignore.stderr = TRUE)

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)

plink <- subset.plink('1KG_EUR/chr' %&&% chr.input,
                      chr.input, ld.lb.input, ld.ub.input, temp.dir)

################################################################
## Take standardized genotype matrix
X <- plink$BED %>% scale() %>% rm.na.zero()
p <- ncol(X)
n.eqtl.size <- ceiling(nrow(X) * sample.ratio)

sim.data <- simulate.gwas(X,
                          eqtl.size = n.eqtl.size,
                          n.med = n.genes,
                          n.causal.med = min(n.genes, n.causal.genes),
                          n.causal.eqtl = n.eqtl.per.gene,
                          n.direct = n.direct.snp,
                          pve.eqtl = h2.eqtl,
                          pve.med = h2.med,
                          pve.direct = h2.unmed)

quick.view <- function(tab, key.var) {
    tab %>% left_join(sim.data$label %>% select(gene, label)) %>%
        arrange_(.dots = key.var) %>% head()
}

################################################################
## Best possible performance (independent mediators)
best.stat <- calc.qtl.stat(sim.data$M.stoch, sim.data$y) %>%
    select(-y.col) %>%
        mutate(best.marg.z = beta/se) %>%
            select(x.col, best.marg.z) %>%
                rename(gene = x.col)

## Best possible performance (multivariate mediators)
best.lm <- lm(scale(sim.data$y) ~ scale(sim.data$M.stoch) - 1)
best.lm.coef <- coefficients(summary(best.lm))

temp <- tibble(gene = 1:n.genes, best.lm.t = best.lm.coef[, 't value'] %>% as.numeric())
best.tab <- best.stat %>% left_join(temp)

best.tab %>% quick.view('desc(abs(best.marg.z))') %>% print()

################################################################
## Imputed TWAS
imp.twas.tab <- run.imputed.twas(X, sim.data)

imp.twas.tab %>% quick.view('desc(abs(imp.lm))') %>% print()

################################################################
## TWAS by Gusev
twas.glmnet.tab <- run.twas.glmnet(X, sim.data)

twas.glmnet.tab %>% quick.view('desc(abs(twas.glmnet.z))') %>% print()

################################################################
## Run TWAS
twas.tab <- run.twas(X, sim.data)

twas.tab %>% quick.view('desc(abs(twas.z))') %>% print()

################################################################
## IVW
ivw.tab <- run.mr(sim.data, X = NULL, is.egger = FALSE) %>%
    mutate(ivw.z = effect/effect.se) %>%
        select(gene, ivw.z)

ivw.tab.rot <- run.mr(sim.data, X = X, is.egger = FALSE) %>%
    mutate(ivw.rot.z = effect/effect.se) %>%
        select(gene, ivw.rot.z)

ivw.tab <-
    ivw.tab %>%
        left_join(ivw.tab.rot)

ivw.tab %>% quick.view('desc(abs(ivw.rot.z))') %>% print()

################################################################
## EGGER
egger.tab <- run.mr(sim.data, X = NULL, is.egger = TRUE) %>%
    mutate(egger.z = effect/effect.se) %>%
        select(gene, egger.z)

egger.tab.rot <- run.mr(sim.data, X = X, is.egger = TRUE) %>%
    mutate(egger.rot.z = effect/effect.se) %>%
        select(gene, egger.rot.z)

egger.tab <-
    egger.tab %>%
        left_join(egger.tab.rot)

egger.tab %>% quick.view('desc(abs(egger.rot.z))') %>% print()

################################################################
## run CaMMEL gene by gene
vb.opt <- list(pi = -1, tau = -5, do.hyper = FALSE,
               tol = 1e-8, vbiter = 3000, do.stdize = TRUE, gammax = 1e4,
               eigen.tol = 1e-2, print.interv = 1000,
               rate = 1e-2, decay = -1e-2, nsample = 10,
               weight = FALSE, do.rescale = FALSE)

z.out <- fit.med.zqtl(effect = sim.data$gwas.beta,
                      effect.se = sim.data$gwas.se,
                      effect.m = sim.data$eqtl.beta,
                      effect.m.se = sim.data$eqtl.se,
                      X.gwas = X,
                      options = vb.opt)

cammel.tab <-
    data.frame(cammel.lodds = z.out$param.med$lodds,
               cammel.z = z.out$param.med$theta / sqrt(z.out$param.med$theta.var)) %>%
                   mutate(gene = 1:n()) %>%
                       arrange(desc(cammel.lodds))

cammel.tab %>% quick.view('desc(cammel.lodds)') %>% print()

out.tab <- sim.data$label %>%
    left_join(cammel.tab) %>%
        left_join(best.stat) %>%
            left_join(imp.twas.tab) %>%
                left_join(twas.tab) %>%
                    left_join(twas.glmnet.tab) %>%
                        left_join(ivw.tab) %>%
                            left_join(egger.tab) %>%
                                arrange(label)

out.tab <- out.tab %>%
    mutate(cammel.lodds = signif(cammel.lodds, digits = 2),
           cammel.z = signif(cammel.z, digits = 2),
           best.marg.z = signif(best.marg.z, digits = 2),
           imp.marg.z = signif(imp.marg.z, digits = 2),
           imp.lm = signif(imp.lm, digits = 2),
           twas.z = signif(twas.z, digits = 2),
           twas.glmnet.z = signif(twas.glmnet.z, digits = 2),
           ivw.z = signif(ivw.z, digits = 2),
           ivw.rot.z = signif(ivw.rot.z, digits = 2),
           egger.z = signif(egger.z, digits = 2),
           egger.rot.z = signif(egger.rot.z, digits = 2),
           ld.idx)

log.msg('Combined all results')

write_tsv(out.tab, path = out.file)

log.msg('Done')
