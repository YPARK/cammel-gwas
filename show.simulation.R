
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
library(PRROC)

sim.files <- list.files('simulation/100_10/', pattern = '*.txt.gz', full.names = TRUE)

.take.num.unmed <- function(s) {
    temp <- basename(s) %>% gsub(replacement = '', pattern = '.txt.gz') %>%
        strsplit(split = '[_]')
    ret <- as.numeric(temp[[1]])
    return(ret[4])
}

sim.num.unmed <- sapply(sim.files, .take.num.unmed) %>%
    unlist(use.names = FALSE)

.read.tsv <- function(j) {
    .sim.types <- paste(c(rep('i',5), rep('d', 17), 'i'), collapse = '')
    ret <- suppressMessages(read_tsv(sim.files[j], col_types = .sim.types))
    if(!('n.unmed' %in% colnames(ret))) {
        ret <- ret %>% mutate(n.unmed = sim.num.unmed[j])
    }
    return(ret)
}

.save.pdf <- function(...) ggsave(..., units = 'in', useDingbats = FALSE, limitsize = FALSE)
.save.png <- function(...) ggsave(..., units = 'in', dpi = 300, limitsize = FALSE)

sim.tab <- lapply(1:length(sim.files), .read.tsv) %>% bind_rows()

sim.tab <- sim.tab %>%
    mutate(CaMMEL = cammel.lodds,
           oTWAS = abs(best.marg.z),
           sTWAS.svd = abs(twas.z),
           sTWAS.glmnet = abs(twas.glmnet.z),
           ivw.2 = abs(ivw.z.2),
           ivw.4 = abs(ivw.z.4),
           ivw.rot.2 = abs(ivw.rot.z.2),
           ivw.rot.4 = abs(ivw.rot.z.4),
           egger.2 = abs(egger.z.2),
           egger.4 = abs(egger.z.4),
           egger.rot.2 = abs(egger.rot.z.2),
           egger.rot.4 = abs(egger.rot.z.4))

################################################################
## Compare AUPRC in general

get.auprc <- function(score, lab) {
    .valid <- !is.na(score)
    if(sum(.valid) < 1) return(NA)
    score <- score[.valid]
    lab <- lab[.valid]
    if(sum(lab == 1) < 1 || sum(lab != 1) < 1) return(NA)
    pr.out <- pr.curve(score[lab == 1], score[lab != 1])
    return(pr.out$auc.davis.goadrich)
}

take.auprc.tab <- function(tab) {
    .ivw <- max(c(get.auprc(tab$ivw.2, tab$label),
                  get.auprc(tab$ivw.4, tab$label),
                  get.auprc(tab$ivw.rot.2, tab$label),
                  get.auprc(tab$ivw.rot.4, tab$label)),
                na.rm = TRUE)

    .egger <- max(c(get.auprc(tab$egger.2, tab$label),
                    get.auprc(tab$egger.4, tab$label),
                    get.auprc(tab$egger.rot.2, tab$label),
                    get.auprc(tab$egger.rot.4, tab$label)),
                  na.rm = TRUE)

    tab %>%
        summarize(CaMMEL = get.auprc(CaMMEL, label),
                  oTWAS = get.auprc(oTWAS, label),
                  sTWAS.svd = get.auprc(sTWAS.svd, label),
                  sTWAS.glmnet = get.auprc(sTWAS.glmnet, label)) %>%
                      mutate(IVW = .ivw, Egger = .egger)
}

methods <- c('oTWAS', 'CaMMEL', 'sTWAS.svd', 'IVW', 'Egger', 'sTWAS.glmnet')
methods.str <- methods %>% gsub(pattern = '[.]', replacement = '-')

method.cols <- c('gray20', 'red', '#4444FF', '#004400', '#009900', '#444499')
method.lty <- c(2, rep(1, 5))
method.shapes <- c(NA, 19, 2, 18, 21, 4)


################################################################
.cond <- c('n.causal.med', 'eqtl.size', 'n.causal.eqtl', 'pve.direct', 'pve.med', 'ld.idx')

auprc.tab <- sim.tab %>%
    group_by_(.dots = .cond) %>%
    do(take.auprc.tab(.)) %>%
    gather_(key_col= 'method', value_col = 'AUPRC', gather_cols = methods) %>%
    group_by_(.dots = c('method', setdiff(.cond, 'ld.idx'))) %>%
    na.omit() %>%
    summarize(AUPRC.mean = mean(AUPRC), AUPRC.se = sd(AUPRC)/sqrt(n()))

df <- auprc.tab %>%
    filter(eqtl.size == 302, pve.direct <= .4, pve.med < .3) %>%
    as.data.frame() %>%
    mutate(pve.direct = (100 * pve.direct) %&&% '% unmediated') %>%
    mutate(n.causal.eqtl = n.causal.eqtl %&&% ' eQTLs / gene') %>%
    mutate(method = factor(method, methods, methods.str))

.aes <- aes(x = 100 * pve.med, y = AUPRC.mean,
            color = method, shape = method, lty = method,
            ymin = (AUPRC.mean - AUPRC.se),
            ymax = (AUPRC.mean + AUPRC.se),
            order = as.numeric(method))

plt <- gg.plot(df, .aes) +
    geom_point() +
    geom_line() +
    geom_linerange(data = df %>% filter(method != 'oTWAS')) +
    scale_linetype_manual(values = method.lty) +
    guides(linetype = guide_legend(nrow = 1)) +
    scale_color_manual(values = method.cols) +
    scale_shape_manual(values = method.shapes) +
    facet_grid(n.causal.eqtl~pve.direct, scales = 'free') +
    ylab('AUPRC') +
    xlab('Proportion of variance explained by mediation (%)') +
    theme(legend.position = 'bottom', legend.title = element_blank())

.save.pdf(filename = 'figure/Fig_sim_auprc.pdf', width = 8, height = 6)
.save.png(filename = 'figure/Fig_sim_auprc.png', width = 8, height = 6)



################################################################
## Compare POWER

get.power <- function(score, lab, fdr.cutoff = 0.01) {

    .valid <- !is.na(score)
    if(sum(.valid) < 2) return(0)
    score <- score[.valid]
    lab <- lab[.valid]
    if(sum(lab == 1) < 1 || sum(lab != 1) < 1) return(0)

    pr.out <- pr.curve(score[lab == 1], score[lab != 1], curve = TRUE)
    .idx <- pr.out$curve[, 2] >= (1 - fdr.cutoff)
    if(sum(.idx) == 0) return(0)
    .temp <- pr.out$curve[.idx, 1]
    return(.temp[1])
}


take.power.tab <- function(tab) {
    .ivw <- max(c(get.power(tab$ivw.2, tab$label),                  
                  get.power(tab$ivw.rot.2, tab$label)),
                na.rm = TRUE)

    .egger <- max(c(get.power(tab$egger.2, tab$label),
                    get.power(tab$egger.rot.2, tab$label)),
                  na.rm = TRUE)

    tab %>%
        summarize(CaMMEL = get.power(CaMMEL, label),
                  oTWAS = get.power(oTWAS, label),
                  sTWAS.svd = get.power(sTWAS.svd, label),
                  sTWAS.glmnet = get.power(sTWAS.glmnet, label)) %>%
                      mutate(IVW = .ivw, Egger = .egger)
}

################################################################
.cond <- c('n.causal.med', 'eqtl.size', 'n.causal.eqtl', 'pve.direct', 'pve.med', 'ld.idx')

.mean <- function(x) mean(x, na.rm = TRUE)
.se <- function(x) sd(x, na.rm = TRUE) / sqrt(pmax(sum(!is.na(x)), 1))

power.tab <- sim.tab %>%
    group_by_(.dots = .cond) %>%
    do(take.power.tab(.)) %>%
    gather_(key_col= 'method', value_col = 'POWER', gather_cols = methods) %>%
    group_by_(.dots = c('method', setdiff(.cond, 'ld.idx'))) %>%
    summarize(POWER.mean = .mean(POWER), POWER.se = .se(POWER))

df <- power.tab %>%
    filter(eqtl.size == 302, pve.direct <= .4, pve.med < .3) %>%
    as.data.frame() %>%
    mutate(pve.direct = (100 * pve.direct) %&&% '% unmediated') %>%
    mutate(n.causal.eqtl = n.causal.eqtl %&&% ' eQTLs / gene') %>%
    mutate(method = factor(method, methods, methods.str))

.aes <- aes(x = 100 * pve.med, y = 100 * POWER.mean,
            color = method, shape = method, lty = method,
            ymin = 100 * (POWER.mean - POWER.se),
            ymax = 100 * (POWER.mean + POWER.se),
            order = as.numeric(method))

plt <- gg.plot(df, .aes) +
    geom_line() +
    geom_linerange(data = df %>% filter(method != 'oTWAS')) +
    geom_point(fill = 'white') +
    scale_linetype_manual(values = method.lty) +
    guides(linetype = guide_legend(nrow = 1)) +
    scale_color_manual(values = method.cols) +
    scale_shape_manual(values = method.shapes) +
    facet_grid(n.causal.eqtl~pve.direct, scales = 'free') +
    ylab('Power (%) at FDR 1%') +
    xlab('Proportion of variance explained by mediation (%)') +
    theme(legend.position = 'bottom', legend.title = element_blank())

.save.pdf(filename = 'figure/Fig_sim_power.pdf', width = 8, height = 6)
.save.png(filename = 'figure/Fig_sim_power.png', width = 8, height = 6)

