
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

sim.files <- list.files('simulation/2', pattern = '*.txt.gz', full.names = TRUE)

.sim.types <- paste(c(rep('i',5), rep('d', 15)), collapse = '')

take.sim.idx <- function(s) {
    temp <- basename(s) %>%
        gsub(pattern='.txt.gz', replacement = '') %>%
            strsplit(split='[_]')
    temp <- temp[[1]]
    as.integer(temp[length(temp)])
}

sim.idx <- sapply(sim.files, take.sim.idx)

.read.tsv <- function(j) {
    suppressMessages(read_tsv(sim.files[j], col_types = .sim.types) %>% mutate(idx = sim.idx[j]))
}

sim.tab <- lapply(1:length(sim.files), .read.tsv) %>% bind_rows()

get.auprc <- function(score, lab) {
    pr.out <- pr.curve(score[lab == 1], score[lab != 1])
    return(pr.out$auc.davis.goadrich)
}

get.power <- function(score, lab, fdr.cutoff = 0.01) {
    pr.out <- pr.curve(score[lab == 1], score[lab == 0], curve = TRUE)
    .idx <- pr.out$curve[, 2] >= (1 - fdr.cutoff)
    if(sum(.idx) == 0) return(0)
    .temp <- pr.out$curve[.idx, 1]
    return(.temp[1])
}

auprc.tab <- function(tab) {
    data.frame(best = get.auprc(tab$best.marg.z, tab$label),
               cammel = get.auprc(tab$cammel.lodds, tab$label),
               twas.glmnet = get.auprc(abs(tab$twas.glmnet.z), tab$label),
               twas.svd = get.auprc(abs(tab$twas.z), tab$label),
               ivw = get.auprc(abs(tab$ivw.z), tab$label),
               ivw.svd = get.auprc(abs(tab$ivw.rot.z), tab$label),
               egger = get.auprc(abs(tab$egger.z), tab$label),
               egger.svd = get.auprc(abs(tab$egger.rot.z), tab$label))
}

power.tab <- function(tab) {
    data.frame(best = get.power(tab$best.marg.z, tab$label),
               cammel = get.power(tab$cammel.lodds, tab$label),
               twas.glmnet = get.power(abs(tab$twas.glmnet.z), tab$label),
               twas.svd = get.power(abs(tab$twas.z), tab$label),
               ivw = get.power(abs(tab$ivw.z), tab$label),
               ivw.svd = get.power(abs(tab$ivw.rot.z), tab$label),
               egger = get.power(abs(tab$egger.z), tab$label),
               egger.svd = get.power(abs(tab$egger.rot.z), tab$label))
}

methods <- c('best', 'cammel', 'twas.glmnet', 'twas.svd', 'ivw', 'ivw.svd', 'egger', 'egger.svd')

temp <- sim.tab %>%
    group_by(n.causal.med, n.causal.eqtl, pve.med, pve.direct, idx) %>%
        do(power.tab(.)) %>%
            gather_(key_col= 'method', value_col = 'power', gather_cols = methods) %>%
                group_by(n.causal.med, n.causal.eqtl, pve.med, pve.direct, method) %>%
                    summarize(power.mu = mean(power), power.se = sd(power)/sqrt(n()))








ggplot(temp %>% filter(n.causal.med == 3),
       aes(x=pve.med, y=power.mu, color = method, lty = method, shape = method)) +
    geom_point() +
    geom_line() +
    facet_grid(n.causal.eqtl ~ pve.direct)
dev.off()


