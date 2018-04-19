################################################################
## construct subnetwork
build.adj <- function(edges) {
    require(Matrix)
    require(methods)

    vertices <- .unlist(edges) %>% unique() %>% sort()

    n <- length(vertices)

    A <- spMatrix(nrow = n, ncol = n,
                  i = match(edges$gene1, vertices),
                  j = match(edges$gene2, vertices),
                  x = rep(1, nrow(edges)))

    A.t <- spMatrix(nrow = n, ncol = n,
                    j = match(edges$gene1, vertices),
                    i = match(edges$gene2, vertices),
                    x = rep(1, nrow(edges)))

    A <- A + A.t
    A[A > 1] <- 1

    D <- spMatrix(nrow = n, ncol = n,
                  i = 1:n, j = 1:n,
                  x = apply(A, 2, sum))

    return(list(A = A, D = D, V = vertices))
}

## remove vertices with iterative degree cutoff
trim.adj <- function(adj, k) {

    A <- adj$A
    D <- adj$D
    V <- adj$V

    rm.genes <- which(diag(D) < k)
    log.msg('Remove %d genes', length(rm.genes))

    while(length(rm.genes) > 0) {

        A <- A[-rm.genes, -rm.genes]
        A[A > 1] <- 1
        V <- V[-rm.genes]

        n <- nrow(A)
        D <- spMatrix(nrow = n, ncol = n,
                      i = 1:n, j = 1:n,
                      x = apply(A, 2, sum))

        rm.genes <- which(diag(D) < k)
        log.msg('Remove %d genes', length(rm.genes))
    }

    list(A = A, D = D, V = V)
}

################################################################
## do spectral clustering by shared neighbor / degree
do.spectral <- function(adj, lambda = 0.01, pve = .9, normalize = TRUE) {

    W <- adj$A
    D <- adj$D
    dd.inv <- diag(D)^(-1/2)

    W <- sweep(sweep(W %*% t(W), 1, dd.inv, `*`), 2, dd.inv, `*`)
    dd <- apply(W, 2, sum)
    dd.inv <- dd^(-1/2)

    n <- length(dd)
    eye <- spMatrix(nrow = n, ncol = n, i = 1:n, j = 1:n, x = rep(1 + lambda, n))
    L <- eye - W
    L.sym <- sweep(sweep(L, 1, dd.inv, `*`), 2, dd.inv, `*`)

    ret <- eigen(as.matrix(L.sym), symmetric = TRUE)

    ## remove negative eigen values
    rm.factors <- which(ret$values <= 1e-8)
    ret$values <- ret$values[-rm.factors]
    ret$vectors <- ret$vectors[, -rm.factors, drop = FALSE]

    ## apply PVE filter
    ll <- ret$values
    k <- max(min(which((cumsum(ll) / sum(ll)) >= pve)), 1)

    log.msg('K = %d features\n', k)

    ret$vectors <- ret$vectors[, 2:(k + 1)]
    ret$values <- ret$values[2:(k = 1)]

    if(normalize) {
        denom <- apply(ret$vectors^2, 1, sum)^(1/2)
        ret$vectors <- sweep(ret$vectors, 1, denom, `/`)
    }

    ret$A <- adj$A
    ret$D <- adj$D
    ret$V <- adj$V
    ret$W <- W

    return(ret)
}

make.cluster <- function(M) {
    require(ClusterR)
    K <- max(ceiling(nrow(M) / 10), 2)
    clust.out <- KMeans_rcpp(M, clusters = K,
                             num_init = 100, max_iters = 100,
                             verbose = FALSE, seed = 1331)
    return(clust.out)
}
