
`%c%` <- function(mat, cols) mat[, cols, drop = FALSE]
`%r%` <- function(mat, rows) mat[rows, , drop = FALSE]
`%&&%` <- function(a,b) paste(a, b, sep = '')
glue <- function(...) paste(..., sep = '')
.unlist <- function(...) unlist(..., use.names = FALSE)
.zeros <- function(n1, n2) matrix(0, n1, n2)

options(stringsAsFactors = FALSE)

.eval <- function(str) eval(parse(text = str))

write.mat <- function(mat, ...) {
    write.table(mat, col.names = FALSE, row.names = FALSE, sep = '\t',
                quote = FALSE, ...)
}

log.msg <- function(...) {
    ss <- as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

.read.mat <- function(...) as.matrix(read.table(...))


fast.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / n.obs
    return(ret)
}

fast.z.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / sqrt(n.obs)
    return(ret)
}

## convert z-score to p-values (two-sided test)
zscore.pvalue <- function(z) {
    2 * pnorm(abs(z), lower.tail = FALSE)
}

## calculate univariate effect sizes and p-values
calc.qtl.stat <- function(xx, yy) {

    require(dplyr)
    require(tidyr)

    .xx <- scale(xx)
    .yy <- scale(yy)

    ## cross-product is much faster than covariance function
    n.obs <- crossprod(!is.na(.xx), !is.na(.yy))
    beta.mat <- crossprod(.xx %>% rm.na.zero(), .yy %>% rm.na.zero()) / n.obs

    log.msg('Computed cross-products')

    ## residual standard deviation
    resid.se.mat <- matrix(NA, ncol(.xx), ncol(.yy))

    for(k in 1:ncol(.yy)) {

        beta.k <- beta.mat[, k]
        yy.k <- .yy[, k]
        err.k <- sweep(sweep(.xx, 2, beta.k, `*`), 1, yy.k, `-`)
        se.k <- apply(err.k, 2, sd, na.rm = TRUE)

        log.msg('Residual on the column %d', k)
        resid.se.mat[, k] <- se.k + 1e-8
    }

    ## organize as consolidated table
    y.cols <- 1:ncol(yy)
    colnames(beta.mat) <- y.cols
    colnames(n.obs) <- y.cols
    colnames(resid.se.mat) <- y.cols

    beta.tab <- beta.mat %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'beta', y.cols)
    
    resid.se.tab <- resid.se.mat %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'resid.se', y.cols)
    
    nobs.tab <- n.obs %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'n', y.cols)
    
    out.tab <- beta.tab %>%
        left_join(nobs.tab) %>%
            left_join(resid.se.tab) %>%
                dplyr::mutate(se = resid.se/sqrt(n)) %>%
                    dplyr::mutate(p.val = zscore.pvalue(beta/se))
    
    out.tab <- out.tab %>%
        mutate(x.col = as.integer(x.col)) %>%
            mutate(y.col = as.integer(y.col))

    return(out.tab)
}

fast.cor <- function(x, y) {
    x.sd <- apply(x, 2, sd, na.rm = TRUE)
    y.sd <- apply(y, 2, sd, na.rm = TRUE)
    ret <- fast.cov(scale(x, scale = FALSE), scale(y, scale = FALSE))
    ret <- sweep(sweep(ret, 1, x.sd, `/`), 2, y.sd, `/`)    
    return(ret)
}

################################################################
## Negative binomial utils
adjust.size.factor <- function(xx) {
    gene.log.mean <- apply(xx, 2, function(x) mean(log(x[x > 0])))
    denom <- as.vector(exp(-gene.log.mean))
    size.factor <- apply(t(xx), 2, function(x) median(x * denom, na.rm = TRUE))
    ret <- apply(xx, 2, function(x) x / size.factor)
    return(ret)
}

stdize.count <- function(xx) {
    .xx <- xx
    .xx[xx <= 0] <- NA
    xx.med <- apply(.xx, 2, median, na.rm = TRUE)
    xx.med <- pmax(xx.med, 1e-4)
    xx.scaled <- sweep(xx, 2, xx.med, `/`)
    ret <- xx.scaled * 50
    return(ret)
}

rm.na.zero <- function(xx) {
    return(replace(xx, is.na(xx), 0))
}

rm.zero <- function(xx) {
    return(replace(xx, xx == 0, NA))
}

trunc <- function(mat, lb = -4, ub = 4) {
    mat[mat > ub] <- ub
    mat[mat < lb] <- lb
    return(mat)
}

################################################################
## Find most correlated (including zero values)
find.cor.idx <- function(Y1, Y0, n.ctrl, p.val.cutoff = 1) {

    colnames(Y1) <- 1:ncol(Y1)
    colnames(Y0) <- 1:ncol(Y0)

    require(dplyr)
    
    y01.stat <- calc.qtl.stat(Y0, Y1) %>%
        dplyr::rename(y0 = x.col, y1 = y.col) %>%
            dplyr::filter(p.val < p.val.cutoff)

    ret <- y01.stat %>% dplyr::group_by(y1) %>%
        dplyr::top_n(n = -n.ctrl, wt = p.val)
    
    return(ret$y0)
}

################################################################
## clean potential genetic signals
subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {
    require(zqtl)
    require(dplyr)
    
    .error <- function(e) {
        log.msg('No QTL here!\n')
        return(NULL)
    }
    
    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {
        
        chr.num <- gsub(pattern = 'chr', replacement = '', chr) %>% as.integer()
        plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr.num, plink.lb, plink.ub, glue(temp.dir, '/plink'))
        system(plink.cmd)
        
        options(stringsAsFactors = FALSE)
        plink <- read.plink(glue(temp.dir, '/plink'))
        colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
        colnames(plink$FAM) <- c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')
        plink$FAM <- plink$FAM %>% mutate(iid = sapply(iid, gsub, pattern = 'GTEX-', replacement = ''))
        
        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 <- 'T'
        }
        
        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 <- 'T'
        }
        return(plink)
    }
    
    plink <- tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error = .error)
    return(plink)
}

################################################################
lift.over.hg38.hg19 <- function(stat.tab, temp.hdr) {

    require(dplyr)
    require(readr)

    dir.create(dirname(temp.hdr), recursive = TRUE)
    temp.in.file <- temp.hdr %&&% '.bed.gz'
    temp.out.file <- temp.hdr %&&% '.lift'
    temp.null.file <- temp.hdr %&&% '.null'

    stat.tab %>%
        dplyr::mutate(snp.loc = as.integer(snp.loc)) %>%
            dplyr::mutate(snp.1 = snp.loc -1) %>%
                dplyr::mutate(prev = snp.loc) %>%
                    dplyr::select(chr, snp.1, snp.loc, rs, prev) %>%
                        write_tsv(path = temp.in.file, col_names = FALSE)

    sys.lift.cmd <- './bin/liftOver' %&&% ' ' %&&% temp.in.file %&&% ' ' %&&% 'hg38ToHg19.over.chain.gz' %&&% ' ' %&&% temp.out.file %&&% ' ' %&&% temp.null.file

    sys.gzip.cmd <- 'gzip -f ' %&&% temp.out.file
    system(sys.lift.cmd %&&% '; ' %&&% sys.gzip.cmd)
    log.msg('System: ' %&&% sys.lift.cmd)
    log.msg('System: ' %&&% sys.gzip.cmd)

    lift.col <- c('chr', 'remove', 'snp.new', 'rs', 'snp.loc')
    lift.tab <- read_tsv(temp.out.file %&&% '.gz', col_names = lift.col) %>%
        dplyr::select(-remove)

    .files <- c(temp.in.file, temp.out.file, temp.out.file %&&% '.gz', temp.null.file)
    unlink(.files)

    ret <- stat.tab %>%
        dplyr::mutate(snp.loc = as.integer(snp.loc)) %>%
            right_join(lift.tab, by = c('chr', 'rs', 'snp.loc')) %>%
                na.omit() %>%
                    dplyr::select(-snp.loc) %>%
                        dplyr::rename(snp.loc = snp.new)

    return(ret)
}

## compare two plink bim files
read.bim <- function(plink.hdr, plink.lb, plink.ub) {
    bim <- read_tsv(plink.hdr %&&% '.bim',
                    col_names = c('chr', 'rs', 'missing', 'snp.loc', 'a1', 'a2'),
                    col_types = 'iciicc')

    bim <- bim %>% filter(snp.loc >= plink.lb, snp.loc <= plink.ub)
    return(bim)
}

match.bim <- function(gwas.bim, qtl.bim) {

    gwas.bim <- gwas.bim %>%
        mutate(gwas.x.pos = 1:n()) %>%
            rename(gwas.plink.a1 = a1,
                   gwas.plink.a2 = a2) %>%
                       select(-missing)

    qtl.bim <- qtl.bim %>%
        mutate(qtl.x.pos = 1:n()) %>%
            rename(qtl.plink.a1 = a1,
                   qtl.plink.a2 = a2,
                   qtl.rs = rs) %>%
                       select(-missing)

    bim.matched <- gwas.bim %>%
        left_join(qtl.bim) %>%
            na.omit()

    return(bim.matched)
}

