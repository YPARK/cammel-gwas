#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
source('Util.R')

library(dplyr)
library(readr)
library(methods)

ld.info.file <- 'ldblocks/EUR/fourier_ls-all.bed'
ld.info.tab <- read_tsv(ld.info.file)

out.dir <- '/broad/hptmp/ypp/cammel/gwas_stat/'
system('mkdir -p ' %&&% out.dir)

## Divide GWAS statistics into independent chunks
write.gwas.chunk <- function(l, gwas.tab, gwas.name) {

    chr.input <- gsub(pattern = 'chr', replacement = '', ld.info.tab[l, 'chr']) %>%
        as.integer()
    ld.lb.input <- ld.info.tab[l, 'start'] %>% as.integer()
    ld.ub.input <- ld.info.tab[l, 'stop'] %>% as.integer()

    out.file <- out.dir %&&% '/' %&&% gwas.name %&&% '_' %&&% l %&&% '.txt.gz'

    if(!file.exists(out.file)) {
        gwas.tab %>%
            dplyr::filter(chr == chr.input, snp.loc >= ld.lb.input, snp.loc <= ld.ub.input) %>%
                dplyr::mutate(gwas.z = gwas.beta / gwas.se) %>%
                    write_tsv(path = out.file)
        log.msg('Written %s', out.file)
        return(1)
    }
    return(0)
}


## load read 1KG variants
bim.cols <- c('chr', 'rs', 'snp.loc')
bim.types <- 'ic_i__'
read.bim <- function(...) read_tsv(..., col_names = bim.cols, col_types = bim.types)

snps.tab <- '1KG_EUR/chr' %&&% 1:22 %&&% '.bim' %>%
    lapply(read.bim) %>%
        bind_rows()

proc.gwas <- function(gwas.file, gwas.name) {

    ptsd.cols <- c('rs', 'gwas.a1', 'gwas.a2', 'gwas.beta', 'gwas.se', 'gwas.p')
    ptsd.types <- 'cccddd_'

    gwas.tab <- read_tsv(gwas.file, col_names = ptsd.cols, col_types = ptsd.types, skip = 1) %>%
        right_join(snps.tab, by = 'rs') %>%
            na.omit() %>%
                mutate(gwas.a1 = toupper(gwas.a1), gwas.a2 = toupper(gwas.a2))

    n.ld <- nrow(ld.info.tab)
    n.written <- sapply(1:n.ld, write.gwas.chunk, gwas.tab = gwas.tab, gwas.name = gwas.name)

    log.msg('Wrote %d files\n', sum(n.written))
    rm(gwas.tab)
    gc()
}

gwas.file <- 'PTSD/military_wave_1.5/ancestry_specific_military/EA/EA_military_7_wave1.51.txt.gz'
gwas.name <- 'ptsd_mil_ea'
proc.gwas(gwas.file, gwas.name)

gwas.file <- 'PTSD/civilian_wave_1.5/ancestry_specific_civilian/EA/EA_civilian_5_wave1.51.txt.gz'
gwas.name <- 'ptsd_civ_ea'
proc.gwas(gwas.file, gwas.name)
