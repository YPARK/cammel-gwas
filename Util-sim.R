
ones <- function(n, m) matrix(1, nrow = n, ncol = m)

zeros <- function(n, m) matrix(0, nrow = n, ncol = m)

.rnorm <- function(n, m) matrix(rnorm(n * m), nrow = n, ncol = m)

match.qtl.plink <- function(qtl.tab, plink.obj) {

    require(dplyr)
    require(tidyr)

    x.bim <- plink.obj$BIM %>%
        dplyr::select(-rs, -missing) %>%
            dplyr::mutate(x.pos = 1:n())

    ret <- qtl.tab %>%
        dplyr::left_join(x.bim, by = c('chr', 'snp.loc')) %>%
            na.omit()

    ret <- ret %>%
        dplyr::filter(((plink.a1 == qtl.a1) & (plink.a2 == qtl.a2)) | ((plink.a1 == qtl.a2) & (plink.a2 == qtl.a1))) %>%
            dplyr::mutate(qtl.z.flip = if_else(qtl.a1 != plink.a1, -qtl.z, qtl.z)) %>%
                dplyr::mutate(qtl.beta.flip = if_else(qtl.a1 != plink.a1, -qtl.beta, qtl.beta))

    ret <- ret %>%
        dplyr::select(-qtl.z, -qtl.beta, -qtl.a1, -qtl.a2) %>%
            rename(qtl.z = qtl.z.flip, qtl.beta = qtl.beta.flip,
                   a1 = plink.a1, a2 = plink.a2)

    return(ret)
}

make.qtl.data <- function(matched.stat) {

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

    .xx <- match(x.pos, qtl.beta$x.pos)
    .mm <- match(med.id, colnames(qtl.beta))
    qtl.beta <- as.matrix(qtl.beta[.xx, .mm])

    .xx <- match(x.pos, qtl.z$x.pos)
    .mm <- match(med.id, colnames(qtl.z))
    qtl.z <- as.matrix(qtl.z[.xx, .mm])

    qtl.se <- qtl.beta / qtl.z
    qtl.se[is.na(qtl.se)] <- mean(qtl.se, na.rm = TRUE)
    qtl.se <- qtl.se + 1e-4

    ret <- list(x.pos = x.pos,
                mediators = med.id,
                qtl.beta = as.matrix(qtl.beta), qtl.se = as.matrix(qtl.se))
    return(ret)
}

check.row.col.order <- function(.tab, rows, cols) {
    .rows <- match(rows, .tab$x.col)
    .tab <- .tab[.rows, -1, drop = FALSE]
    .cols <- match(cols, colnames(.tab))
    .tab <- .tab[, .cols, drop = FALSE]
}

simulate.gwas <- function(X,
                          eqtl.size,
                          n.med,
                          n.causal.med,
                          n.causal.eqtl,
                          n.direct,
                          pve.eqtl,
                          pve.med,
                          pve.direct){

    n.ind <- nrow(X)
    p <- ncol(X)

    ## 0. sample residuals
    var.resid <- pmax(1 - pve.med - pve.direct, 1e-8)
    log.msg('Residual Var = %.2e', var.resid)
    y <- .rnorm(n.ind, 1) * sqrt(var.resid)

    print(var(y))

    ## uniformly distributed genes
    gene.loc <- sample(p, n.med) %>% sort()

    ## 1. sample causal eQTL SNPs
    M <- matrix(NA, nrow = n.ind, ncol = n.med)
    causal.eqtl.tab <- NULL
    sd.eqtl <- sqrt(pve.eqtl/n.causal.eqtl)

    for(k in 1:n.med) {

        dist.k <- abs(1:p - gene.loc[k])
        weight.k <- exp(-dist.k / n.med)

        alpha.k <- .rnorm(n.causal.eqtl, 1) * sd.eqtl
        causal.snp.k <- sample(p, n.causal.eqtl, prob = weight.k)
        .df <- data.frame(snp = causal.snp.k, gene = k, alpha = alpha.k)
        causal.eqtl.tab <- bind_rows(causal.eqtl.tab, .df)
        M[, k] <- X %c% causal.snp.k %*% alpha.k
    }

    M.err <- .rnorm(n.ind, n.med) * sqrt(1 - pve.eqtl)

    ## 2. sample causal mediation genes
    med.var <- pve.med / n.causal.med
    log.msg('Mediation Var = %.2e', med.var)

    causal.med.tab <- data.frame(gene = sample(n.med, n.causal.med),
                                 beta = as.numeric(.rnorm(n.causal.med, 1) * sqrt(med.var)))

    causal.genes <- unique(causal.med.tab$gene)

    M.stoch <- M + M.err

    y <- y + (M.stoch %c% causal.med.tab$gene) %*% matrix(causal.med.tab$beta, n.causal.med, 1)

    causal.med.tab <-
        causal.med.tab %>%
            left_join(causal.eqtl.tab)

    ## 3. sample direct effect
    pleiotropic.tab <- causal.eqtl.tab %>%
        filter(!(gene %in% causal.genes)) %>%
            filter(!(snp %in% causal.med.tab$snp))

    ## Select non-causal pleiotropic genes, trying to achieve InSIDE
    wrong.genes <- pleiotropic.tab %>%
        sample_n(n.direct)

    if(n.direct > 2) {
        gam.mat <- .rnorm(n.direct, 1000)
        abs.cor <- abs(apply(gam.mat, 2, cor, y = wrong.genes$alpha, use = 'pairwise.complete.obs'))
        gam <- gam.mat %c% which.min(abs.cor)
    } else {
        gam <- rnorm(n.direct)
    }

    wrong.genes <- wrong.genes %>%
        mutate(gamma = gam)

    inside.cor <- cor(wrong.genes$alpha, wrong.genes$gamma) %>% as.numeric()

    eta.direct <- (X %c% as.integer(wrong.genes$snp)) %*% matrix(wrong.genes$gamma, ncol = 1)
    sig.direct <- sd(eta.direct)
    eta.direct <- eta.direct / sig.direct * sqrt(pve.direct)

    y <- y + eta.direct

    ## calculate summary statistics
    eqtl.size <- min(eqtl.size, nrow(M.stoch))
    eqtl.idx <- sample(nrow(M.stoch), eqtl.size)
    X.eqtl <- X %r% eqtl.idx
    M.eqtl <- M.stoch %r% eqtl.idx

    eqtl.M <- calc.qtl.stat(X.eqtl, M.eqtl)

    eqtl.beta <- eqtl.M %>% select(beta, x.col, y.col) %>%
        spread(key = y.col, value = beta) %>%
            check.row.col.order(rows = 1:p, cols = 1:n.med)

    eqtl.se <- eqtl.M %>% select(se, x.col, y.col) %>%
        spread(key = y.col, value = se) %>%
            check.row.col.order(rows = 1:p, cols = 1:n.med)

    gwas <- calc.qtl.stat(X, y)

    gwas.beta <- gwas %>% select(beta, x.col, y.col) %>%
        spread(key = y.col, value = beta) %>%
            check.row.col.order(rows = 1:p, cols = 1)

    gwas.se <- gwas %>% select(se, x.col, y.col) %>%
        spread(key = y.col, value = se) %>%
            check.row.col.order(rows = 1:p, cols = 1)

    inside.cor <- signif(inside.cor, 2)

    label <- rep(0, n.med)
    label[unique(wrong.genes$gene)] <- -1
    label[unique(causal.med.tab$gene)] <- 1
    label.tab <- data.frame(gene = 1:n.med, label = label) %>%
        mutate(eqtl.size,
               n.causal.med,
               n.causal.eqtl,
               pve.eqtl,
               pve.med,
               pve.direct,
               inside.cor) %>%
                   as.data.frame()

    ret <- list(med = causal.med.tab,
                wrong.med = wrong.genes,
                eqtl = causal.eqtl.tab,
                eqtl.beta = as.matrix(eqtl.beta),
                eqtl.se = as.matrix(eqtl.se),
                gwas.beta = as.matrix(gwas.beta),
                gwas.se = as.matrix(gwas.se),
                M = M,
                M.stoch = M.stoch,
                y = y,
                M.eqtl = M.eqtl,
                X.eqtl = X.eqtl,
                label = label.tab)

    log.msg('Generated simulation data')
    return(ret)
}

################################################################
## run gene by gene TWAS by Mancuso et al. 2017
run.twas <- function(X, sim.data, do.stdize = TRUE) {

    if(do.stdize) {
        gwas.z <- (sim.data$gwas.beta / sim.data$gwas.se) %>% as.matrix()
        eqtl.Z <- (sim.data$eqtl.beta / sim.data$eqtl.se) %>% as.matrix()
    } else {
        gwas.z <- (sim.data$gwas.beta) %>% as.matrix()
        eqtl.Z <- (sim.data$eqtl.beta) %>% as.matrix()
    }

    n.med <- ncol(eqtl.Z)

    svd.out <- take.ld.svd(X, options = list(do.stdize = TRUE, eigen.tol = 1e-2))
    W.t <- sweep(svd.out$V.t, 1, svd.out$D, `/`)
    ## V.t <- sweep(svd.out$V.t, 1, svd.out$D, `*`)
    ## LD.inv <- t(W.t) %*% W.t
    ## LD <- t(V.t) %*% V.t

    take.twas <- function(g) {

        qtl.z <- eqtl.Z %c% g

        ## gwas.mult <- t(W.t) %*% (W.t %*% gwas.z)
        ## qtl.mult <- t(W.t) %*% (W.t %*% qtl.z)
        ## mult <- sum(qtl.mult * gwas.mult) / sqrt(sum(qtl.mult^2) + 1e-8)

        ## This is too memory-intensive
        ## qtl.z.poly <- LD.inv %*% qtl.z
        ## num <- t(qtl.z.poly) %*% gwas.z
        ## denom <- t(qtl.z.poly) %*% LD %*% qtl.z.poly

        num <- t(W.t %*% gwas.z) %*% (W.t %*% qtl.z)
        denom <- t(W.t %*% qtl.z) %*% (W.t %*% qtl.z)

        ret <- signif(as.numeric(num/sqrt(denom + 1e-16)), 4)
        log.msg('TWAS finished [%d / %d] : %.2f', g, n.med, ret)

        data.frame(twas.z = ret)
    }

    twas.result <- bind_rows(lapply(1:n.med, take.twas))

    twas.tab <- twas.result %>%
        mutate(gene = 1:n()) %>%
            arrange(desc(abs(twas.z)))

    log.msg('TWAS Finished\n\n')
    return(twas.tab)
}

## Run TWAS of Gusev using glmnet
run.twas.glmnet <- function(X, sim.data, do.stdize = TRUE) {

    n.med <- ncol(sim.data$M.eqtl)    
    M.train <- sim.data$M.eqtl %>% scale()
    X.train <- sim.data$X.eqtl %>% scale()

    X.test <- X %>% scale() %>% rm.na.zero()
    yy.test <- sim.data$y %>% scale()

    ## sTWAS = dot(z.GWAS, z.true.eQTL) / sqrt(z.true.eQTL' * R * z.true.eQTL + 1e-16)

    if(do.stdize) {
        gwas.z <- (sim.data$gwas.beta / sim.data$gwas.se) %>% as.matrix()
    } else {
        gwas.z <- (sim.data$gwas.beta) %>% as.matrix()
    }

    svd.out <- take.ld.svd(X, options = list(do.stdize = TRUE, eigen.tol = 1e-2))
    Vd.t <- sweep(svd.out$V.t, 1, svd.out$D, `*`)

    .twas.glmnet.fun <- function(g) {
        .yy <- M.train %c% g
        effect <- run.glmnet(.yy, X.train, alpha = 0.5)

        if(sum(abs(effect)) < 1e-8) return(data.frame(twas.glmnet.z = 0)) # no eQTL found

        .snps <- which(abs(effect) > 0)
        .xx <- X.train %c% .snps
        .df <- data.frame(.yy, .xx)
        colnames(.df) <- c('y', as.character(.snps))
        .lm <- lm(y ~ . - 1, data = .df)
        .eqtl.z <- coefficients(summary(.lm))
        .snps <- rownames(.eqtl.z) %>% gsub(pattern = '\`', replacement = '') %>%
            as.integer()
        .eqtl.z <- .eqtl.z %c% 't value' %>% as.matrix()

        num <- t(gwas.z %r% .snps) %*% .eqtl.z %>% as.numeric()
        eta.g <- Vd.t %c% .snps %*% .eqtl.z
        denom <- sqrt(sum(eta.g^2) + 1e-16)
        ret <- data.frame(twas.glmnet.z = signif(num/denom, 4))
        return(ret)
    }

    n.med <- ncol(M.train)

    twas.result <- bind_rows(lapply(1:n.med, .twas.glmnet.fun))

    twas.tab <- twas.result %>%
        mutate(gene = 1:n()) %>%
            arrange(desc(abs(twas.glmnet.z)))

    log.msg('TWAS Finished\n\n')
    return(twas.tab)
}

run.glmnet <- function(y, x, alpha = 1){
    require(glmnet)
    require(methods)
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])
    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    ret <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    pnz <- mean(abs(ret) > 0)
    log.msg('non-zeros = %.2f', pnz)
    return(ret)
}

run.imputed.twas <- function(X, sim.data) {

    n.med <- ncol(sim.data$M.eqtl)    
    M.train <- sim.data$M.eqtl %>% scale()
    X.train <- sim.data$X.eqtl %>% scale()

    X.test <- X %>% scale() %>% rm.na.zero()
    yy.test <- sim.data$y %>% scale()

    take.imp <- function(g) {
        effect <- run.glmnet(M.train %c% g, X.train, alpha = 0.5)
        m.imp <- X.test %*% effect %>% as.numeric()
        ret <- data.frame(g = m.imp)
        colnames(ret) <- g
        return(ret)
    }

    M.imp <- lapply(1:n.med, take.imp) %>% bind_cols()

    ## independent mediators
    imp.stat <- calc.qtl.stat(M.imp, yy.test) %>%
        select(-y.col) %>%
            mutate(imp.marg.z = beta/se) %>%
                select(x.col, imp.marg.z) %>%
                    rename(gene = x.col)
    
    ## Best possible performance (another glmnet
    imp.lm <- run.glmnet(yy.test, M.imp %>% as.matrix())
    imp.lm.coef <- imp.lm %>% unlist(use.names = FALSE) %>% as.numeric()
    
    temp <- tibble(gene = 1:n.med, imp.lm = imp.lm.coef)
    imp.stat <- imp.stat %>% left_join(temp)
    
    return(imp.stat)
}

################################################################
## Run MR-Egger (optinally with PC weighting)
run.mr <- function(sim.data, X = NULL, is.egger = TRUE, weak.cutoff = 0) {
    require(MendelianRandomization)
    require(zqtl)

    eqtl.beta <- sim.data$eqtl.beta
    eqtl.se <- sim.data$eqtl.se
    gwas.beta <- sim.data$gwas.beta
    gwas.se <- sim.data$gwas.se

    if(!is.null(X)) {
        ## Rotate to make them orthogonal
        svd.out <- take.ld.svd(X, options = list(do.stdize = TRUE, eigen.tol = 1e-2))
        W.t <- sweep(svd.out$V.t, 1, svd.out$D, `/`)

        eqtl.beta <- W.t %*% (eqtl.beta / eqtl.se)
        eqtl.se <- matrix(1, nrow(eqtl.beta), ncol(eqtl.beta))

        gwas.beta <- W.t %*% (gwas.beta / gwas.se)
        gwas.se <- matrix(1, nrow(gwas.beta), ncol(gwas.beta))
    }

    build.mr.obj <- function(j) {
        z.j <- eqtl.beta[, j] / eqtl.se[, j]
        valid <- which(abs(z.j) > weak.cutoff)
        if(length(valid) < 3) return(NULL)
        mr.obj <- mr_input(bx = eqtl.beta[valid, j],
                           bxse = eqtl.se[valid, j],
                           by = gwas.beta[valid, 1],
                           byse = gwas.se[valid, 1])
    }

    take.egger <- function(j) {
        mr.obj <- build.mr.obj(j)
        if(is.null(mr.obj)) return(data.frame(effect = 0, effect.se = 1e-4))

        egger.out <- mr_egger(mr.obj)
        data.frame(effect = egger.out@Estimate,
                   effect.se = egger.out@StdError.Est)
    }

    take.ivw <- function(j) {
        mr.obj <- build.mr.obj(j)
        if(is.null(mr.obj)) return(data.frame(effect = 0, effect.se = 1e-4))

        ivw.out <- mr_ivw(mr.obj)
        data.frame(effect = ivw.out@Estimate,
                   effect.se = ivw.out@StdError)
    }

    fun <- take.egger
    if(!is.egger) fun <- take.ivw

    n.med <- eqtl.beta %>% ncol()

    mr.result <- lapply(1:n.med, fun) %>%
        bind_rows() %>%
            mutate(gene = 1:n.med) %>%
                arrange(desc(abs(effect/effect.se)))

    return(mr.result)
}
