
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
## match allele directions in statistics
## input
## 1. plink
## 2. gwas.tab
## 3. qtl.tab
match.allele <- function(gwas.tab, plink.obj, qtl.tab) {
    require(dplyr)

    x.bim <- plink.obj$BIM %>%
        dplyr::select(-rs, -missing) %>%
            dplyr::mutate(x.pos = 1:n())

    ret <- qtl.tab %>%
        dplyr::left_join(x.bim, by = c('chr', 'snp.loc')) %>%
            dplyr::left_join(gwas.tab, by = c('chr', 'snp.loc')) %>%
                na.omit()

    ret <- ret %>%
        dplyr::filter(((plink.a1 == qtl.a1) & (plink.a2 == qtl.a2)) | ((plink.a1 == qtl.a2) & (plink.a2 == qtl.a1))) %>%
            dplyr::mutate(gwas.z.flip = if_else(qtl.a1 != plink.a1, -gwas.z, gwas.z)) %>%
                dplyr::mutate(gwas.beta.flip = if_else(qtl.a1 != plink.a1, -gwas.beta, gwas.beta)) %>%
                    dplyr::mutate(qtl.z.flip = if_else(qtl.a1 != plink.a1, -qtl.z, qtl.z)) %>%
                        dplyr::mutate(qtl.beta.flip = if_else(qtl.a1 != plink.a1, -qtl.beta, qtl.beta))
    
    ret <- ret %>%
        dplyr::select(-gwas.z, -gwas.beta, -qtl.z, -qtl.beta,
                      -qtl.a1, -qtl.a2, -gwas.a1, -gwas.a2) %>%
                          rename(gwas.z = gwas.z.flip, gwas.beta = gwas.beta.flip,
                                 qtl.z = qtl.z.flip, qtl.beta = qtl.beta.flip,
                                 a1 = plink.a1, a2 = plink.a2)

    return(ret)
}



## match evething to plink.gwas
match.plink <- function(plink.gwas, plink.qtl) {

    if(is.null(plink.gwas)) return(NULL)
    if(is.null(plink.qtl)) return(NULL)

    ret.gwas <- plink.gwas
    ret.qtl <- plink.qtl

    gwas.bim <- plink.gwas$BIM %>%
        mutate(gwas.x.pos = 1:n()) %>%
            rename(gwas.plink.a1 = plink.a1,
                   gwas.plink.a2 = plink.a2) %>%
                       select(-missing)

    qtl.bim <- plink.qtl$BIM %>%
        mutate(qtl.x.pos = 1:n()) %>%
            rename(qtl.plink.a1 = plink.a1,
                   qtl.plink.a2 = plink.a2,
                   qtl.rs = rs) %>%
                       select(-missing)

    bim.matched <- gwas.bim %>%
        left_join(qtl.bim) %>%
            na.omit()

    if(nrow(bim.matched) < 1) return(NULL)

    bim.matched <- bim.matched %>%
        dplyr::filter(((gwas.plink.a1 == qtl.plink.a1) & (gwas.plink.a2 == qtl.plink.a2)) |
                          ((gwas.plink.a2 == qtl.plink.a1) & (gwas.plink.a1 == qtl.plink.a2))) %>%
                              arrange(chr, snp.loc)

    if(nrow(bim.matched) < 1) return(NULL)

    ret.gwas$BIM <- ret.gwas$BIM[bim.matched$gwas.x.pos, , drop = FALSE]
    ret.gwas$BED <- ret.gwas$BED[ , bim.matched$gwas.x.pos, drop = FALSE]

    ret.qtl$BIM <- ret.qtl$BIM[bim.matched$qtl.x.pos, , drop = FALSE]
    ret.qtl$BED <- ret.qtl$BED[ , bim.matched$qtl.x.pos, drop = FALSE]

    flip.tab <- ret.gwas$BIM %>% mutate(gwas.x.pos = 1:n()) %>%
        left_join(ret.qtl$BIM %>% mutate(qtl.x.pos = 1:n()),
                  by = c('chr', 'snp.loc'),
                  suffix = c('.gwas', '.qtl')) %>%                      
                      filter(plink.a1.gwas != plink.a1.qtl)

    ret.qtl$BIM[flip.tab$qtl.x.pos, ] <- ret.gwas$BIM[flip.tab$gwas.x.pos, ]

    flip.bed <- ret.qtl$BED[, flip.tab$qtl.x.pos]
    zero.idx <- flip.bed <= 0.5
    two.idx <- flip.bed >= 1.5
    flip.bed[two.idx] <- 0
    flip.bed[zero.idx] <- 2
    ret.qtl$BED[, flip.tab$qtl.x.pos] <- flip.bed

    return(list(gwas = ret.gwas, qtl = ret.qtl))
}

################################################################
## generate matrices from matched statistics
make.zqtl.data <- function(matched.stat) {

    require(dplyr)
    require(tidyr)

    if(nrow(matched.stat) < 1) {
        return(NULL)
    }

    qtl.beta <- matched.stat %>% dplyr::select(med.id, x.pos, qtl.beta) %>%
        dplyr::group_by(med.id, x.pos) %>%
            dplyr::slice(which.max(abs(qtl.beta))) %>%
                tidyr::spread(key = med.id, value = qtl.beta)

    qtl.z <- matched.stat %>% dplyr::select(med.id, x.pos, qtl.z) %>%
        dplyr::group_by(med.id, x.pos) %>%
            dplyr::slice(which.max(abs(qtl.z))) %>%
                tidyr::spread(key = med.id, value = qtl.z)

    med.id <- colnames(qtl.beta)[-1]
    x.pos <- qtl.beta$x.pos

    .temp <- matched.stat %>% dplyr::select(x.pos, gwas.beta, gwas.se) %>%
        group_by(x.pos) %>% dplyr::slice(which.max(abs(gwas.beta / gwas.se)))

    gwas.beta <- qtl.beta %>%
        dplyr::select(x.pos) %>%
            left_join(.temp %>% dplyr::select(x.pos, gwas.beta))

    gwas.se <- qtl.beta %>%
        dplyr::select(x.pos) %>%
            left_join(.temp %>% dplyr::select(x.pos, gwas.se))

    .xx <- match(x.pos, qtl.beta$x.pos)
    .mm <- match(med.id, colnames(qtl.beta)) 
    qtl.beta <- as.matrix(qtl.beta[.xx, .mm])

    .xx <- match(x.pos, qtl.z$x.pos)
    .mm <- match(med.id, colnames(qtl.z)) 
    qtl.z <- as.matrix(qtl.z[.xx, .mm])

    qtl.se <- qtl.beta / qtl.z
    qtl.se[is.na(qtl.se)] <- mean(qtl.se, na.rm = TRUE)
    qtl.se <- qtl.se + 1e-4

    .xx <- match(x.pos, gwas.beta$x.pos)
    gwas.beta <- as.matrix(gwas.beta[.xx, 'gwas.beta'])

    .xx <- match(x.pos, gwas.se$x.pos)
    gwas.se <- as.matrix(gwas.se[.xx, 'gwas.se'])

    gwas.se[is.na(gwas.se)] <- mean(gwas.se, na.rm = TRUE)

    ret <- list(x.pos = x.pos,
                mediators = med.id,
                qtl.beta = as.matrix(qtl.beta), qtl.se = as.matrix(qtl.se),
                gwas.beta = as.matrix(gwas.beta), gwas.se = as.matrix(gwas.se))
    return(ret)
}

################################################################
## genrate theoretical null data

make.zqtl.null <- function(X, se.mat, eig.tol, beta.mat = NULL, is.indep = TRUE, stdize = TRUE) {
    require(zqtl)

    .rnorm <- function(nrow, ncol) matrix(rnorm(nrow * ncol), nrow, ncol)

    n.traits <- ncol(se.mat)

    ## Estimate LD matrix
    svd.out <- take.ld.svd(X, options = list(do.stdize = stdize, eigen.tol = eig.tol))

    if(n.traits < 2 || is.indep) {

        ## R = V' D2 V
        ## (a) sample theta.tilde ~ N(0, D^2) <=> D * N(0, I)
        ## (b) sample beta.hat ~ se * V * theta.tilde 
        d <- nrow(svd.out$D)
        theta.tilde <- sweep(.rnorm(d, n.traits), 1, svd.out$D, `*`)
        beta.null.mat <- (t(svd.out$V.t) %*% theta.tilde) * se.mat

    } else {

        Vd <- sweep(t(svd.out$V.t), 2, svd.out$D, `*`) # R = Vd * Vd'
        d <- ncol(Vd)

        ## Estimate correlation between z-scores
        z.mat <- beta.mat / se.mat
        z.cov <- cor(z.mat, use = 'pairwise.complete.obs') # scale of covariance is too large
        L <- chol(z.cov, pivot = FALSE) # C = L' * L

        z.null <- Vd %*% .rnorm(d, n.traits) %*% L
        beta.null.mat <- z.null * se.mat 

    }

    return(beta.null.mat)
}

################################################################
## Estimate variance
estimate.variance <- function(z.out, z.data, X, nn) {

    V.t <- z.out$Vt
    Y <- z.out$Y
    M <- z.out$M
    S.inv <- z.out$S.inv.y
    S <- 1/S.inv
    S.inv.m <- z.out$S.inv.m
    D2 <- z.out$D2
    D <- sqrt(z.out$D2)
    kk <- length(D)

    Vd <- sweep(t(V.t), 2, D, `*`)
    W <- sweep(t(V.t), 2, D, `/`)

    ## Estimate total variance using Shi et al.
    eta.gwas <- t(W) %*% z.data$gwas.beta
    gg <- sum(eta.gwas^2)
    var.shi <- (gg * nn - kk) / (nn - kk)

    ## Estimated true qtl effect
    ## theta.hat        ~ S R inv(S) (aa * bb)
    ## inv(S) theta.hat ~ R inv(S) (aa * bb)
    ## -- the mediated component
    aa <- sweep(W %*% sweep(M, 1, D, `/`), 1, S, `*`)
    bb <- z.out$param.mediated$theta
    ab <- aa %*% bb
    ## -- direct effect
    theta.dir <- z.out$param.direct$theta

    xx.std <- scale(X) %>% rm.na.zero()
    n.ref <- nrow(xx.std)

    ## residual variance
    r.hat <- z.out$resid$theta
    rr <- sweep(W %*% (t(W) %*% r.hat), 1, S, `*`)
    eta.r <- t(Vd) %*% rr
    var.resid <- sum(eta.r^2)

    ## Estimate variance using reference genotype matrix
    eta.tot <- xx.std %*% (theta.dir + ab)
    eta.dir <- xx.std %*% (theta.dir)
    eta.med.tot <- xx.std %*% (ab)

    var.ref.tot <- var(eta.tot, na.rm = TRUE) + var.resid
    var.ref.dir <- var(eta.dir, na.rm = TRUE)
    var.ref.med <- var(eta.med.tot, na.rm = TRUE)

    ## Estimate variance using covariance
    eta.tot <- t(Vd) %*% (theta.dir + ab)
    var.tot <- sum(eta.tot^2)
    eta.dir <- t(Vd) %*% (theta.dir)
    var.dir <- sum(eta.dir^2)
    eta.ab <- t(Vd) %*% ab
    var.med <- sum(eta.ab^2)

    ## each mediation effect
    var.each <- function(k) {
        ab.k <- (aa %c% k) %*% (bb %r% k)
        eta.ab.k <- t(Vd) %*% ab.k
        var.k <- sum(eta.ab.k^2)
        return(var.k)
    }

    n.med <- nrow(bb)
    var.med.vec <- sapply(1:n.med, var.each)

    return(list(shi = signif(var.shi, digits = 4),
                resid = signif(var.resid, 4),
                ref.tot = signif(var.ref.tot, 4),
                ref.dir = signif(var.ref.dir, 4),
                ref.med = signif(var.ref.med, 4),
                tot = signif(var.tot, 4),
                dir = signif(var.dir, 4),
                med.tot = signif(var.med, 4),
                med = signif(var.med.vec, 4),
                k = kk))
}

################################################################
melt.effect <- function(effect.obj, .rnames, .cnames) {
    require(reshape2)
    melt.list <- lapply(seq_along(effect.obj),
                        function(mat.list, val.names, i) {
                            mat <- signif(mat.list[[i]], 4)
                            val <- val.names[[i]]
                            rownames(mat) <- .rnames
                            colnames(mat) <- .cnames
                            melt(mat, value.name = val) },
                        mat.list = effect.obj, val.names = names(effect.obj))

    ret <- Reduce(function(...) left_join(..., by = c('Var1', 'Var2')), melt.list) %>%
        mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
    return(ret)
}

################################################################
row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D <- proxy::dist(mat, method = function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] <- 0
    h.out <- hclust(D)
    o.out <- cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

gg.plot <- function(...) {
    require(ggplot2)
    ggplot(...) + theme_bw() + theme(plot.background = element_blank(),
                                     panel.background = element_blank(),
                                     strip.background = element_blank(),
                                     legend.background = element_blank())
}

order.pair <- function(pair.tab) {

    require(tidyr)
    require(dplyr)
    M <- pair.tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro <- row.order(M %>% dplyr::select(-row) %>% as.matrix())
    log.msg('Sort the rows: %d', length(ro))
    co <- row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    log.msg('Sort the columns: %d', length(co))
    cc <- colnames(M)[-1]
    rr <- M[, 1]

    list(rows = rr[ro], cols = cc[co], M = M)
}
