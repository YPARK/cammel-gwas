#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-Figure.R')

library(dplyr)
library(tidyr)
library(readr)
library(methods)
library(ggplot2)
library(ggrepel)
library(circlize)


################################################################
draw.circos <- function(med.tab, out.file, fsr.cutoff = 1e-2) {

    pdf(out.file, useDingbats = FALSE, width = 7, height = 7)

    circos.initializeWithIdeogram(plotType = c('labels', 'axis'), species = 'hg19',
                                  chromosome.index = 'chr' %&&% 1:22)

    fun.gwas <- function(x, y) {
        .gwas <- -log10(5e-8)
        .sig <- which(y > .gwas)
        .null <- which(y <= .gwas)
        circos.lines(x = c(min(x),max(x)), y = c(.gwas, .gwas), col = 'red', lty = 2)
        circos.points(x[.null], y[.null], pch = 19, cex = .1, col = 'gray60')
        circos.points(x[.sig], y[.sig], pch = 19, cex = .5, col = 'gray40')
    }
    
    circos.par('track.height' = .05)
    circos.track(factors = 'chr' %&&% med.tab$chr,
                 x = med.tab$gene.loc,
                 bg.border = 'gray80',
                 bg.lwd = 0.1,
                 y = pmin(-log10(med.tab$gwas.p.by.gwas), 10),
                 panel.fun = fun.gwas)

    ## show gene density
    circos.par('track.height' = .1)
    circos.trackHist(factors = 'chr' %&&% med.tab$chr,
                     x = med.tab$gene.loc, bin.size = 5e6,
                     bg.border = 'gray80',
                     bg.lwd = 0.1,
                     col = 'gray80')

    circos.yaxis('left', labels.cex = .3, sector.index = 'chr1', col = 'orange', labels.col = 'orange')

    ## show mediation FSR
    fun.lfsr <- function(x, y) {
        .fsr <- -log10(fsr.cutoff)
        .sig <- which(y > .fsr)
        .null <- which(y <= .fsr)
        circos.lines(x = c(min(x),max(x)), y = c(.fsr,.fsr), col = 'red', lty = 2)
        circos.points(x[.null], y[.null], pch = 16, cex = .1, col = 'gray40')
        circos.points(x[.sig], y[.sig], pch = 16, cex = .5, col = 'gray20')
    }

    circos.par('track.height' = .1)
    circos.track(factors = 'chr' %&&% med.tab$chr,
                 x = med.tab$gene.loc,
                 y = -log10(med.tab$lfsr),
                 bg.border = 'gray90',
                 bg.col = 'gray90',
                 panel.fun = fun.lfsr)

    circos.yaxis('left', labels.cex = .3, sector.index = 'chr1', col = 'orange', labels.col = 'orange')

    ## show mediation effect size
    fun.theta <- function(x, y) {
        .up <- which(y > 1e-3)
        .dn <- which(y < -1e-3)
        .null <- which(abs(y) < 1e-3)
        circos.points(x[.null], y[.null], pch = 19, cex = .1, col = 'gray40')
        circos.points(x[.up], y[.up], pch = 19, cex = .5, col = 'red')
        circos.points(x[.dn], y[.dn], pch = 19, cex = .5, col = 'blue')
    }

    circos.par('track.height' = .15)
    circos.track(factors = 'chr' %&&% med.tab$chr,
                 x = med.tab$gene.loc,
                 y = med.tab$theta,
                 bg.col = 'gray90',
                 bg.border = 'gray90',
                 panel.fun = fun.theta)

    circos.yaxis('left', labels.cex = .3, sector.index = 'chr1', col = 'orange', labels.col = 'orange')

    ## show top 3 gene names within each LD block
    med.lab.tab <-
        med.tab %>%
            filter(lfsr <= fsr.cutoff) %>%
                group_by(chr, ld.lb, ld.ub) %>%
                    top_n(3, -var.mediated) %>%
                        as.data.frame() %>%
                            mutate(gene.col = if_else(theta > 0, 'red', 'blue')) %>%
                                select(chr, tss, tes, hgnc, gene.col, lfsr, var.mediated) %>%
                                    unique() %>%
                                        arrange(chr) %>%
                                            mutate(chr = 'chr' %&&% chr)

    ## at most top 70 genes
    med.lab.tab <- med.lab.tab %>% arrange(desc(var.mediated)) %>%
        head(70)

    circos.par('track.height' = .2)
    circos.genomicLabels(med.lab.tab,
                         labels.column = 4,
                         side = 'inside',
                         padding = .5,
                         cex = .4,
                         col = med.lab.tab$gene.col,
                         line_col = 'gray50')

    circos.clear()
    dev.off()
}

read.med.tab <- function(med.file, fsr.cutoff = 1e-2) {
    med.tab <-
        read_tsv(med.file) %>%
            mutate(gene.loc = if_else(strand == '+', tss, tes)) %>%
                mutate(lfsr = pmax(lfsr, fsr.cutoff/10)) %>%
                    group_by(hgnc) %>%
                        slice(which.min(lfsr)) %>%
                            as.data.frame()
}

med.tab <- read.med.tab('mediation/gene_ad_mayo_gammax-4_eigen-2.txt.gz', 5e-2)
out.file <- 'figure/gene/ad/Fig_global_mayo.pdf'
draw.circos(med.tab, out.file, 5e-2)

med.tab <- read.med.tab('mediation/gene_ad_rosmap_gammax-4_eigen-2.txt.gz', 5e-2)
out.file <- 'figure/gene/ad/Fig_global_rosmap.pdf'
draw.circos(med.tab, out.file, 5e-2)

med.tab <- read.med.tab('mediation/gene_ad_geuvadis_gammax-4_eigen-2.txt.gz', 5e-4)
out.file <- 'figure/gene/ad/Fig_global_geuvadis.pdf'
draw.circos(med.tab, out.file, 1e-3)
