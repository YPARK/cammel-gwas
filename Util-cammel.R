run.cammel <- function(zqtl.data, xx.gwas, xx.med, opt) {

    if(is.null(zqtl.data)) return(NULL)

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()

    z.out <- fit.med.zqtl(zqtl.data$gwas.beta, zqtl.data$gwas.se,
                          zqtl.data$qtl.beta, zqtl.data$qtl.se,
                          X.gwas = xx.gwas.mat,
                          X.med = xx.med.mat,                          
                          options = opt)
    return(z.out)
}

run.cammel.null <- function(zqtl.data, xx.gwas, xx.med, n.null, opt) {

    xx.gwas.mat <- xx.gwas %c% zqtl.data$x.pos %>% as.matrix()
    xx.med.mat <- xx.med %c% zqtl.data$x.pos %>% as.matrix()
    ret <- NULL
    for(r in 1:n.null) {

        gwas.null <- make.zqtl.null(xx.gwas.mat, zqtl.data$gwas.se, eig.tol = eig.tol)
        qtl.null <- make.zqtl.null(xx.med.mat, zqtl.data$qtl.se, eig.tol = eig.tol)

        z.out <- fit.med.zqtl(gwas.null, zqtl.data$gwas.se,
                              qtl.null, zqtl.data$qtl.se,
                              X.gwas = xx.gwas.mat, X.med = xx.med.mat,                              
                              options = opt)

        null.stat <- melt.effect(z.out$param.mediated, zqtl.data$mediators, r) %>%
            mutate(theta.var = sqrt(theta.var)) %>% rename(theta.se = theta.var) %>%
                mutate(log.var = as.numeric(signif(log(z.out$var.decomp$var.med.each), 2)))

        print(null.stat %>% filter(lodds > 0) %>% as.data.frame())
        ret <- bind_rows(ret, null.stat)
        log.msg('\nnull round = %d / %d\n', r, n.null)
    }
    ret <- ret %>% rename(med.id = Var1, null = Var2)
    return(ret)
}

get.var.tab <- function(var.decomp, mediators) {

    if(is.null(var.decomp)) return(NULL)

    ret <- data.frame(med.id = mediators,
                      var.mediated = signif(var.decomp$var.med.each, 2),
                      var.mediated.tot = signif(var.decomp$var.med.mean, 2),
                      var.mediated.tot.se = signif(sqrt(var.decomp$var.med.var), 2),
                      var.direct.tot = signif(var.decomp$var.direct.mean, 2),
                      var.direct.tot.se = signif(sqrt(var.decomp$var.direct.var), 2))
    return(ret)
}

## gene-level QTL and GWAS stat summary
get.summary.tab <- function(.gwas.tab, .qtl.tab, .plink.obj) {
    if(is.null(.gwas.tab)) return(NULL)

    if('rs' %in% colnames(.qtl.tab)){
        .qtl.tab <- .qtl.tab %>% select(-rs)
    }

    .temp <- .gwas.tab %>%
        match.allele(plink.obj = .plink.obj, qtl.tab = .qtl.tab)

    .temp.gwas <- .temp %>% group_by(med.id) %>%
        slice(which.min(gwas.p)) %>%
            select(med.id, rs, snp.loc, gwas.p, gwas.beta, gwas.se, qtl.z)

    .temp.qtl <- .temp %>% group_by(med.id) %>%
        slice(which.max(qtl.z)) %>%
            select(med.id, rs, snp.loc, gwas.p, gwas.beta, gwas.se, qtl.z)

    ret <- .temp.gwas %>%
        left_join(.temp.qtl, by = 'med.id', suffix = c('.by.gwas', '.by.qtl'))

    ret <- ret %>% dplyr::select_(.dots = sort(names(ret)))

    return(ret)
}

get.effect.tab <- function(z.out, z.data, gwas.tab, qtl.tab, data.name, plink.obj = NULL) {
    if(is.null(z.out)) return(NULL)
    z.effect <-
        melt.effect(z.out$param.mediated, z.data$mediators, data.name) %>%
            rename(med.id = Var1, gwas = Var2) %>%
                left_join(get.var.tab(z.out$var.decomp, z.data$mediators))

    if(!is.null(plink.obj)) {
        z.effect <- z.effect %>%
            left_join(get.summary.tab(gwas.tab, qtl.tab, plink.obj))
    }
    return(z.effect)
}
