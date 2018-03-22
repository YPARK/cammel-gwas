#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
source('Util.R')
source('Util-Figure.R')

library(zqtl)
library(dplyr)
library(readr)
library(methods)
library(ggplot2)
library(ggrepel)

################################################################
## Functions to make the code more readable
read.igap.data <- function(argv) {

    if(length(argv) < 8) {
        return(NULL)
    }

    ld.idx <- as.integer(argv[1])          # e.g., ld.idx = 736
    qtl.file <- argv[2]                    # e.g., qtl.file = 'cis-eqtl/rosmap/736_qtl.txt.gz'
    qtl.sample.size <- as.integer(argv[3]) # e.g., qtl.sample.size = 500
    geno.dir <- argv[4]                    # e.g., geno.dir = 'ROSMAP_GENO' # (eQTL genotype matrix)
    gammax.input <- as.numeric(argv[5])    # e.g., gammax.input = 1e4
    eig.tol <- as.numeric(argv[6])         # e.g., eig.tol = 1e-2
    med.file <- argv[7]                    # e.g., med.file = 'mediation/gene_ad_rosmap_gammax-4_eigen-2.txt.gz'
    out.hdr <- argv[8]                     # e.g., out.hdr = 'temp.cammel'

    dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

    ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
    ld.info.tab <- read_tsv(ld.info.file)
    chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[ld.idx, 'chr']) %>%
        as.integer()
    ld.lb.input <- ld.info.tab[ld.idx, 'start'] %>% as.integer()
    ld.ub.input <- ld.info.tab[ld.idx, 'stop'] %>% as.integer()

    med.tab <-
        read_tsv(med.file) %>%
            filter(ld.lb >= ld.lb.input, ld.ub <= ld.ub.input, chr == chr.input) %>%
                as.data.frame()

    gwas.dir <- './gwas_stat/'
    igap.gwas.file <- gwas.dir %&&% '/igap_' %&&% ld.idx %&&% '.txt.gz'

    .files <- c(igap.gwas.file)
    if(!all(sapply(.files, file.exists))) {
        log.msg('Missing input files: %s\n', paste(.files, collapse = ', '))
        return(NULL)
    }

    if(!file.exists(qtl.file)) {
        log.msg('QTL file does not exist: %s\n', qtl.file)
        return(NULL)
    }

    genes <- med.tab %>%
        filter(lfsr < 1e-2) %>%
            select(hgnc) %>%
                unlist(use.names = FALSE) %>%
                    paste(collapse = '_')
    
    out.hdr <- out.hdr %&&% '_' %&&% chr.input %&&% '_' %&&% ld.idx %&&% '_' %&&% genes

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
        return(NULL)
    }

    plink.gwas <-  plink.matched$gwas
    plink.eqtl <-  plink.matched$qtl

    ## Read QTL statistics and measure basic statistics
    ## chr rs snp.loc med.id qtl.a1 qtl.a2 qtl.beta qtl.z
    ## i   c  i       c      c      c      d        d
    qtl.tab <- read_tsv(qtl.file, col_types = 'icicccdd')

    if(nrow(qtl.tab) == 0) {
        log.msg('Empty QTL file\n')
        return(NULL)
    }

    qtl.tab <- qtl.tab %>% dplyr::select(-rs)

    med.qtl.stat <- qtl.tab %>%
        group_by(med.id) %>%
            summarize(n.qtl.4 = sum(abs(qtl.z) > 4))

    igap.gwas.tab <- read_tsv(igap.gwas.file)

    igap.matched <- igap.gwas.tab %>%
        match.allele(plink.obj = plink.eqtl, qtl.tab = qtl.tab)

    igap.data <- igap.matched %>%
        make.zqtl.data()

    gc()

    ret <- list(plink.gwas = plink.gwas,
                plink.eqtl = plink.eqtl,
                zqtl.data = igap.data,
                gammax = gammax.input,
                eig.tol = eig.tol,
                n.eqtl = qtl.sample.size,
                mediation = med.tab,
                gwas.matched = igap.matched,
                out.hdr = out.hdr)

    system('rm -r ' %&&% temp.dir)
    log.msg('Successfully read the data!\n')
    return(ret)
}

################################################################
## Combine z-scores
build.zscore.df <- function(zqtl.out, x.bim) {

    z.obs <- t(zqtl.out$Vt) %*% zqtl.out$Y %>% as.numeric()

    z.unmed.mat <- zqtl.out$resid.Z

    z.unmed <- apply(z.unmed.mat, 1, mean, na.rm = TRUE) %>%
        matrix(ncol=1) %*%
            zqtl.out$param.unmediated$theta %>%
                as.numeric()

    z.unmed.sd <- apply(z.unmed.mat, 1, sd, na.rm = TRUE) %>%
        matrix(ncol=1) %*%
            zqtl.out$param.unmediated$theta %>%
                as.numeric()

    z.med.hat <- t(zqtl.out$Vt) %*% zqtl.out$M %*% zqtl.out$param.mediated$theta %>%
        as.numeric()

    ret <- tibble(snp.loc = x.bim$snp.loc,
                  obs = z.obs,
                  unmed = z.unmed,
                  unmed.sd = z.unmed.sd,
                  med = z.med.hat)

    ret <- ret %>% mutate(resid = obs - unmed)
}

rotate.zscore <- function(zscore.df, z.out) {

    Vt <- z.out$Vt
    dd <- sqrt(z.out$D2)

    ret <- tibble(factor = 1:nrow(Wt),
                  lambda = dd %>% as.numeric(),
                  obs = Vt %*% zscore.df$obs %>% as.numeric(),
                  unmed = Vt %*% zscore.df$unmed %>% as.numeric(),
                  resid = Vt %*% zscore.df$resid %>% as.numeric(),
                  med = Vt %*% zscore.df$med %>% as.numeric())
}

################################################################
## show eQTL edges
build.eqtl.tab <- function(mediation, lodds.cutoff = 0) {

    med.info <- mediation %>%
        select(med.id, tss, tes, strand, lodds) %>%
            as.data.frame()

    eqtl.z.tab <-
        (igap.data$zqtl.data$qtl.beta / igap.data$zqtl.data$qtl.se) %>%
            as.data.frame() %>% mutate(x = 1:n()) %>%
                gather(key = 'med.id', value = 'z', -x) %>%
                    left_join(x.bim) %>%
                        left_join(med.info) %>%
                            filter(lodds > lodds.cutoff) %>%
                                select(tss, tes, snp.loc, strand, z) %>%
                                    mutate(gene.loc = if_else(strand == '+', tss, tes)) %>%
                                        as.data.frame()

    return(eqtl.z.tab)
}

build.gene.tab <- function(med.tab) {
    ret <- med.tab %>%
        mutate(gene.loc = if_else(strand == '+', tss, tes)) %>%
            mutate(gene.end = if_else(strand == '+', tes, tss)) %>%
                arrange(desc(lodds)) %>%
                    select(hgnc, gene.loc, gene.end, strand)

    ret <- ret %>%
        mutate(y.pos = 1:n()) %>%
            mutate(hj = if_else(strand == '+', 1, 0))

    return(ret)
}

l10p <- function(z) -log10(pmax(2 * pnorm(abs(z), lower.tail = FALSE), 1e-10))

################################################################
## convert z-score matrix to p-value tibble
convert.zscore.pval <- function(z.df) {

    z.unmed <- z.df$unmed
    z.unmed.sd <- z.df$unmed.sd

    p.lb <- (z.unmed + 2*z.unmed.sd) %>% l10p()
    p.ub <- (z.unmed - 2*z.unmed.sd) %>% l10p()

    ret <- tibble(snp.loc = z.df$snp.loc,
                  obs = z.df$obs %>% l10p(),
                  unmed = z.unmed %>% l10p(),
                  unmed.max = pmax(p.lb, p.ub),
                  med = z.df$med %>% l10p(),
                  resid = z.df$resid %>% l10p())
}


################################################################
## To make x-axis scale
build.x.scale <- function(plink.bim.tab, x.pos, mediation.tab = NULL) {

    kb.lab <- function(x) format(round(x/1e3), big.mark = ',')

    .bim.tab <- plink.bim.tab %r% x.pos %>%
        mutate(x = 1:n())

    if(!is.null(mediation.tab)) {
        x.min <- min(c(.bim.tab$snp.loc, mediation.tab$tss))
        x.max <- max(c(.bim.tab$snp.loc, mediation.tab$tes))
    } else {
        x.min <- min(.bim.tab$snp.loc)
        x.max <- max(.bim.tab$snp.loc)
    }

    x.scale <- scale_x_continuous(limits = c(x.min, x.max),
                                  labels = kb.lab)

    x.scale.top <- scale_x_continuous(limits = c(x.min, x.max),
                                      labels = kb.lab,
                                      position = 'top')

    x.grid <- seq(x.min, x.max, length = nrow(.bim.tab))

    ret <- list(scale = x.scale,
                scale.top = x.scale.top,
                bim.tab = .bim.tab,
                grid = x.grid)
}

################################################################
##
build.ld.cov <- function(xx.gwas.mat, xx.bim, xx.grid, ...) {

    ld.pairs <- take.ld.pairs(xx.gwas.mat, ...) %>%
        as.data.frame()

    .temp <- tibble(x = 1:length(xx.grid), x.grid = xx.grid)

    x.coord <- xx.bim %>% filter(x %in% unique(ld.pairs$x)) %>%
        select(snp.loc, x) %>%
            left_join(.temp)

    n.snps <- nrow(x.coord)
    x.range <- max(x.coord$x.grid) - min(x.coord$x.grid)

    x.min <- min(xx.grid)

    ld.tab <-
        ld.pairs %>% left_join(x.coord) %>%
            mutate(x.pos = x.min + (x.pos - 0.5) * x.range/n.snps) %>%
                mutate(g = x %&&% '_' %&&% y)

    ret <- list(ld.cov = ld.tab, grid2loc = x.coord)
}


################################################################
## Run the method again to estimate effects
run.cammel <- function(data.list) {

    igap.sample.size <- 74056
    qtl.sample.size <- data.list$n.eqtl
    gammax.input <- data.list$gammax
    eig.tol <- data.list$eig.tol

    plink.eqtl <- data.list$plink.eqtl
    plink.gwas <- data.list$plink.gwas
    zqtl.data <- data.list$zqtl.data

    vb.opt <- list(pi.ub = -1, pi.lb = -5, tau = -5, do.hyper = TRUE, tol = 1e-8, gammax = gammax.input,
                   vbiter = 3000, do.stdize = TRUE, eigen.tol = eig.tol,
                   rate = 1e-2, nsample = 10, print.interv = 1500,
                   weight = FALSE, do.rescale = FALSE)

    xx.gwas <- plink.gwas$BED
    xx.med <- plink.eqtl$BED
    gwas.sample.size <- igap.sample.size

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()

    z.out <- fit.med.zqtl(zqtl.data$gwas.beta, zqtl.data$gwas.se,
                          zqtl.data$qtl.beta, zqtl.data$qtl.se,
                          X.gwas = xx.gwas.mat,
                          X.med = xx.med.mat,
                          n = gwas.sample.size,
                          n.med = qtl.sample.size, options = vb.opt)
}

.grid.arr <- function(x, ...) grid.arrange(grobs = x, padding = unit(0, 'in'), ...)
.save.pdf <- function(...) ggsave(..., units = 'in', useDingbats = FALSE, limitsize = FALSE)
.save.png <- function(...) ggsave(..., units = 'in', dpi = 300, limitsize = FALSE)

################################################################
## argv <- c('736',
##           'cis-eqtl/rosmap/736_qtl.txt.gz',
##           '500',
##           'ROSMAP_GENO',
##           '1e4',
##           '1e-2',
##           'mediation/gene_ad_rosmap_gammax-4_eigen-2.txt.gz',
##           'temp.cammel')

igap.data <- read.igap.data(argv)

if(is.null(igap.data)) { q() }

med.out.file <- igap.data$out.hdr %&&% '_mediation.pdf'

if(file.exists(med.out.file)) { q() }

x.sc <- build.x.scale(igap.data$plink.gwas$BIM, igap.data$zqtl.data$x.pos, igap.data$mediation)
x.sc.nomed <- build.x.scale(igap.data$plink.gwas$BIM, igap.data$zqtl.data$x.pos)
x.bim <- x.sc$bim.tab
z.out <- run.cammel(igap.data)

################################################################
## Fig 1. show z-scores on the genomic axis
zscore.df <- build.zscore.df(z.out, x.bim)

.aes <- aes(x = snp.loc,
            y = unmed,
            ymin = unmed - 2*unmed.sd,
            ymax = unmed + 2*unmed.sd)

.df <-
    zscore.df %>%
        gather(key = 'data', value = 'z', obs, unmed, resid) %>%
            mutate(unmed.sd = if_else(data == 'unmed', unmed.sd, 0)) %>%
                as.data.frame()

## truncate between -4 and 4 for visualization
.df <- .df %>% mutate(z = pmax(z, -4)) %>% mutate(z = pmin(z, 4))

.df$data <- factor(.df$data,
                   c('obs', 'unmed', 'resid'),
                   c('observed GWAS', 'unmediated effects',
                     'residuals w/o the unmediated'))

plt.z1 <- gg.plot(.df) +
    geom_ribbon(data = .df %>% filter(unmed.sd > 0),
                aes(x = snp.loc, ymin = z - unmed.sd, ymax = z + unmed.sd),
                fill = 'gray60') +
                    geom_point(aes(x = snp.loc, y = z, color = data), size = .3) +
                        scale_color_manual(values = c('gray20', '#006600', '#000066'), guide = FALSE) +
facet_wrap(~data, nrow = 1) +
    x.sc.nomed$scale + ylab('z-score') + xlab('genomic locations (kb)') +
        theme(strip.text.x = element_text(size = 8),
              axis.text.x = element_text(size = 5))

plt.z2 <- gg.plot(.df) +
    geom_point(aes(x = med, y = z, color = data), size = .3) +
        geom_smooth(aes(x = med, y = z), method = 'lm', se = FALSE, size = .5, color = 'red') +
            scale_color_manual(values = c('gray20', '#006600', '#000066'), guide = FALSE) +
facet_wrap(~data, nrow = 1) +
    ylab('z-score') + xlab('estimated mediation z-score') +
        theme(strip.text.x = element_blank(),
              axis.text = element_text(size = 8))

grid.out.zscore <-
    list(plt.z1, plt.z2) %>% match.widths() %>%
        .grid.arr(ncol = 1)

out.file <- igap.data$out.hdr %&&% '_zscore.pdf'
.save.pdf(filename = out.file, plot = grid.out.zscore, width = 6, height = 5)

out.file <- igap.data$out.hdr %&&% '_zscore.png'
.save.png(filename = out.file, plot = grid.out.zscore, width = 6, height = 5)

################################################################
## Fig 2. show rotated z-scores
rot.df <- rotate.zscore(zscore.df, z.out)

.df <-
    rot.df %>%
    gather(key = 'data', value = 'z', obs, unmed, resid) %>%
    as.data.frame()

.df$data <- factor(.df$data,
                   c('obs', 'unmed', 'resid'),
                   c('observed GWAS', 'unmediated effects',
                     'residuals w/o the unmediated'))

plt.rot <- 
    gg.plot(.df, aes(x = med, y = z, color = data)) +
    geom_point(aes(size = lambda)) +
    geom_smooth(method = 'lm', se = FALSE, size = .5, color = 'red') +
    scale_size_continuous(range = c(0, 2), guide = FALSE) +
    scale_color_manual(values = c('gray20', '#006600', '#000066'), guide = FALSE)

plt.rot <- plt.rot +
    facet_wrap(~data, nrow = 1, scales = 'free') +
    ylab('rotated z-score') + xlab('rotated mediation z-score') +
    theme(strip.text.x = element_blank(),
          axis.text = element_text(size = 8))

out.file <- igap.data$out.hdr %&&% '_rotated.pdf'
.save.pdf(filename = out.file, plot = plt.rot, width = 6, height = 2.5)

out.file <- igap.data$out.hdr %&&% '_rotated.png'
.save.png(filename = out.file, plot = plt.rot, width = 6, height = 2.5)

################################################################
## Fig 3. show p-values with LD
pval.df <- convert.zscore.pval(zscore.df)

ld.list <- build.ld.cov(igap.data$plink.gwas$BED %c% igap.data$zqtl.data$x.pos,
                        x.bim,
                        x.sc.nomed$grid,
                        stdize = TRUE,
                        cutoff = 0.25)

ld.tab <- ld.list$ld.cov
grid2loc <- ld.list$grid2loc

plt.bridge <-
    gg.plot(grid2loc, aes(x = snp.loc, xend = x.grid, y = 1, yend = 0)) +
        geom_segment(color = 'gray', alpha = 0.5, size = .5) +
            x.sc.nomed$scale.top +
                theme_void()

plt.ld <-
    gg.plot(ld.tab) +
        geom_polygon(aes(x = x.pos, y = y.pos, fill = abs(cov), group = g)) +
            scale_fill_continuous(low = 'white', high = 'black', guide = FALSE) +
                theme_void()

plt.med <-
    gg.plot(pval.df) +
        geom_point(aes(x = snp.loc, y = med), color = 'gray20', size = .5) +
            x.sc.nomed$scale + xlab('genomic location (kb)') +
                ylab('mediated (-log10 P)')

grob.med <-
    list(plt.med, plt.bridge, plt.ld) %>%
        match.widths()

grid.out.ld <- grob.med %>% .grid.arr(heights = c(2, .5, 1.5))

ww <- (max(x.sc$grid) - min(x.sc$grid)) / 1e6 * 2 + .5

out.file <- igap.data$out.hdr %&&% '_ld.png'
.save.png(filename = out.file, plot = grid.out.ld, width = ww, height = 6)


################################################################

scale.col <- scale_color_gradient2(low = '#222288', high = '#882222', mid = 'gray', guide = FALSE)

plt.gwas <-
    gg.plot(igap.data$gwas.matched, aes(x = snp.loc, y = pmin(-log10(gwas.p), 10))) +
        geom_point(size = .5, color = 'gray20') + x.sc$scale +
            theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
                ylab('total GWAS')

med.tab <- igap.data$mediation %>%
    mutate(gene.loc = if_else(strand == '+', tss, tes))

plt.med.fsr <-
    gg.plot(med.tab, aes(x = gene.loc, y = -log10(lfsr))) +
    geom_hline(yintercept = 2, lty = 2, color = 'orange') +
    geom_segment(stat = 'identity', aes(xend = gene.loc, color = theta), yend = 0, size = 2) +
    x.sc$scale + scale.col +
    xlab('genomic location (kb)') + ylab('-log10 FSR (mediation)')

med.tab.sig <- med.tab %>% filter(lfsr < .5)

plt.med.fsr <- 
    plt.med.fsr +
    geom_text_repel(data = med.tab.sig,
                    aes(label = hgnc %&&% ' (' %&&% signif(lfsr,2) %&&% ')'),
                    nudge_y = .2, nudge_x = .1, size = 3, color = 'gray20', segment.color = 'gray')

plt.med <-
    gg.plot(pval.df) +
        geom_point(aes(x = snp.loc, y = med), color = 'gray20', size = .5) +
            x.sc$scale.top + scale.col + xlab('genomic location (kb)') +
                ylab('mediated (-log10 P)')

eqtl.z.tab <- build.eqtl.tab(igap.data$mediation) %>%
    mutate(z = pmax(z, -2)) %>%
        mutate(z = pmin(z, 2))

z.color.scale <- scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                                       limits = c(-2, 2), guide = FALSE)

plt.eqtl <- gg.plot(eqtl.z.tab %>% filter(abs(z) > 1)) +
    theme_void() +
        geom_segment(aes(x = gene.loc, xend = snp.loc, color = z), y = 1, yend = 0) +
            z.color.scale +
                x.sc$scale

gene.aes <- aes(x = gene.loc, xend = gene.end, y = y.pos, yend = y.pos,
                label = hgnc, hjust = hj)
gene.arrow <- arrow(length = unit(.5, 'lines'))
gene.df <- build.gene.tab(med.tab)
gene.df.sig <- build.gene.tab(med.tab %>% filter(lfsr < .5))
gene.df.long <- gene.df %>% filter(abs(gene.end - gene.loc) > 1e5)

plt.gene <-
    gg.plot(gene.df, gene.aes) +
    geom_text_repel(data = gene.df.sig, size = 3) +
    geom_point(size = 1, pch = 15, color = '#005500')

plt.gene <- plt.gene +
    x.sc$scale + theme_void() + scale_y_continuous(expand = c(.5, 0))

plt.gene <- plt.gene +
    geom_segment(data = gene.df.long, arrow = gene.arrow, size = 1, color = '#005500')

plt.gene <- plt.gene +
    geom_point(data = gene.df.sig, color = 'red', pch = 21, size = 1.5)

grid.out.med <-
    list(plt.med.fsr, plt.gene, plt.eqtl, plt.med + scale_y_reverse(), plt.gwas + scale_y_reverse()) %>%
        match.widths() %>%
            .grid.arr(heights = c(3, 1, .7, 2, 1.8), ncol = 1)

ww <- (max(x.sc$grid) - min(x.sc$grid)) / 1e6 * 2 + .5

out.file <- igap.data$out.hdr %&&% '_mediation.png'
.save.png(filename = out.file, plot = grid.out.med, width = ww, height = 7)

out.file <- igap.data$out.hdr %&&% '_mediation.pdf'
.save.pdf(filename = out.file, plot = grid.out.med, width = ww, height = 7)
