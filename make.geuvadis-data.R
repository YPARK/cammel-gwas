#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(dplyr)
library(readr)
library(methods)
source('Util.R')

expr.file <- 'GEUVADIS_EXPR/GD660.GeneQuantCount.txt.gz'
covar.file <- '1KG_EUR/1kg_sample_info.txt.gz'

expr.dat <- read_tsv(expr.file)

genes.info <- expr.dat[, 1:4] %>%
    dplyr::select(-Gene_Symbol) %>%
        rename(ensg = TargetID, chr = Chr, gene.loc = Coord)

Y <- expr.dat[, -(1:4)] %>% t()

.temp <- apply(Y < 1, 2, mean, na.rm = TRUE)
.rm <- which(is.na(.temp) | .temp > .5)

if(length(.rm) > 0) {
    Y <- Y[, -.rm, drop = FALSE]
    genes.info <- genes.info[ -.rm, , drop = FALSE]
}

sample.info <- tibble(Sample = colnames(expr.dat[, -(1:4)])) %>%
    dplyr::mutate(Sample = substr(Sample, start = 1, stop = 7)) %>%
        mutate(idx = 1:n()) %>%
            group_by(Sample) %>%
                summarize(n = n(), pos = list(idx))

Y.mean <- matrix(NA, nrow(sample.info), ncol(Y))
for(r in 1:nrow(sample.info)) {
    r.idx <- unlist(sample.info[r, ]$pos) 
    Y.mean[r, ] <- apply(Y[r.idx, , drop = FALSE], 2, mean, na.rm = TRUE)    
}

Y.mean <- Y.mean %>% adjust.size.factor()
Y.out <- t(Y.mean) %>% as.data.frame()
colnames(Y.out) <- sample.info$Sample

dir.create('cis-eqtl/geuvadis_data', recursive = TRUE)

write_tsv(genes.info, path = 'cis-eqtl/geuvadis_data/genes.txt.gz')
write_tsv(sample.info %>% select(Sample, n), path = 'cis-eqtl/geuvadis_data/samples.txt.gz')
write_tsv(Y.out, path = 'cis-eqtl/geuvadis_data/expression.txt.gz')

