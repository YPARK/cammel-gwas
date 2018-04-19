#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

if(length(argv) != 2) q()

library(dplyr)
library(readr)
library(goseq)
library(tidyr)
library(ggrepel)
source('Util.R')
source('Util-clustering.R')

med.file <- argv[1] # e.g., 'mediation/gene_igap_rosmap_gammax-4_eigen-2.txt.gz'
out.hdr <- argv[2]  # e.g., ''

dir.create(dirname(out.hdr), recursive = TRUE)

################################################################
cutoff <- 0.01
max.sz.cutoff <- 1000

################################################################
med.tab <- read_tsv(med.file) %>%
    dplyr::group_by(med.id) %>%
    dplyr::slice(which.min(lfsr)) %>%
    as.data.frame()

net.file <- 'pathways/BIOGRID-ORGANISM-Homo_sapiens-3.4.159.tab2.txt.gz'
net.tab <- read_tsv(net.file) %>%
    dplyr::filter(`Experimental System Type` == 'physical') %c%
    c(8, 9) %>% unique()
colnames(net.tab) <- c('gene1', 'gene2')

.save.pdf <- function(...) ggsave(..., units = 'in', useDingbats = FALSE, limitsize = FALSE)
.save.png <- function(...) ggsave(..., units = 'in', dpi = 300, limitsize = FALSE)

################################################################
## 0. Takek GO slim

library(biomaRt)
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
                  path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

go.slim <- getBM(filters = 'ensembl_gene_id',
                 values = med.tab$med.id,
                 attributes = 'goslim_goa_accession',
                 mart = grch37)

go <- getgo(med.tab$med.id, 'hg19', 'ensGene')
go.slim.cat <- lapply(go, function(x) x[x %in% go.slim[, 1]])

################################################################
## 1. Take significant genes.  These genes will be seed genes.
sig.hgnc <- med.tab %>% dplyr::filter(lfsr < cutoff) %>%
    dplyr::select(hgnc) %>% unique() %>% .unlist()

sig.ensg <- med.tab %>% dplyr::filter(lfsr < cutoff) %>%
    dplyr::select(med.id) %>% unique() %>% .unlist()

################################################################
## 2. Construct subnetworks from the seed, looking at immediate neighbors.
net.genes <- net.tab %>% dplyr::filter(gene1 %in% sig.hgnc | gene2 %in% sig.hgnc) %>%
    dplyr::select(gene1, gene2) %>% .unlist() %>% unique()

net.genes <- c(net.genes, sig.hgnc) %>% unique() %>% sort()

net.edges <- net.tab %>% dplyr::filter(gene1 %in% net.genes, gene2 %in% net.genes) %>%
    dplyr::filter(gene1 != gene2) %>%
    unique()

################################################################
## 3. Build edges then apply iterative degree cutoff >= 2
trimmed <- net.edges %>% build.adj() %>% trim.adj(k = 2)

sig.genes <- med.tab %>% dplyr::filter(hgnc %in% trimmed$V | lfsr < cutoff) %>%
    dplyr::select(med.id) %>% unique() %>% .unlist()

genes.tot <- med.tab %>% dplyr::select(med.id) %>% unique() %>% .unlist()
genes.vector <- as.integer(genes.tot %in% sig.genes)
names(genes.vector) <- genes.tot

pwf <- nullp(genes.vector, 'hg19', 'ensGene', plot.fit = FALSE)

go.enrich.tab <- goseq(pwf, genome = 'hg19', id = 'ensGene', gene2cat = go.slim.cat) %>%
    dplyr::mutate(p.val = pmin(over_represented_pvalue, 1)) %>%
    dplyr::select(category, p.val, ontology, term, numDEInCat, numInCat) %>%
    mutate(q.val = p.adjust(p.val, 'fdr')) %>%
    arrange(q.val)

################################################################
## 4. Show enrichment of all the genes connected to the significant
## genes

build.gene.go.tab <- function(show.genes) {
    genes2go <- getgo(show.genes$med.id, genome = 'hg19', id = 'ensGene')
    go.list <- genes2go[!is.na(names(genes2go))]
    n.genes <- length(go.list)
    genes <- names(go.list)
    take.go.j <- function(j) data.frame(med.id = genes[j], category = go.list[[j]])
    ret <- lapply(1:n.genes, take.go.j) %>% bind_rows() %>%
        left_join(show.genes)
    return(ret)
}

order.genes.terms <- function(tab) {
    temp <- tab %>%
        dplyr::select(term, hgnc, lodds) %>%
            dplyr::rename(col = term, row = hgnc, weight = lodds) %>%
                order.pair()

    co <- match(temp$cols, names(temp$M))
    M <- temp$M[, co]
    ro <- order(apply(abs(M) > 0, 1, function(x) median(which(x))), decreasing = TRUE)
    terms <- temp$cols %>% as.character()
    genes <- temp$M$row[ro] %>% as.character()
    return(list(terms = terms, genes = genes))
}



show.genes <- med.tab %>% dplyr::filter(hgnc %in% net.genes) %>%
    dplyr::select(med.id, hgnc, lodds, theta, theta.var, lfsr)    

gene.go.tab <- build.gene.go.tab(show.genes)

write_tsv(gene.go.tab, path = out.hdr %&&% '_tab.txt.gz')

go.cat <- gene.go.tab$category %>% unique()

sig.genes <- med.tab %>% dplyr::filter(lfsr < cutoff) %>%
    dplyr::select(med.id, hgnc, lodds, theta, theta.var, lfsr)    

gene.go.tab.sig <- build.gene.go.tab(sig.genes) # GO categories containing significant genes

go.tab.show <- gene.go.tab %>% left_join(go.enrich.tab) %>%
    na.omit() %>%
    as.data.frame() %>% 
    dplyr::filter(q.val < 1e-1, numInCat <= max.sz.cutoff) %>%
    na.omit() %>%
    as.data.frame()

orders <- go.tab.show %>% 
    order.genes.terms()

tab <- go.tab.show %>%
    mutate(term = factor(term, orders$terms),
           hgnc = factor(hgnc, orders$genes))

tab.sig <- tab %>% dplyr::filter(lfsr < cutoff) %>%
    mutate(z.show = pmin(pmax(theta/sqrt(theta.var), -4), 4))

tab.sig.annot <- tab.sig %>% dplyr::group_by(hgnc) %>%
    dplyr::slice(which.min(p.val))

rr <- tab$category %>% unique() %>% length()
cc <- tab$hgnc %>% unique() %>% length()

ww <- cc / 300 + sapply(tab$term, function(x) nchar(as.character(x))) %>% max() / 10
hh <- rr / 9 + 1

p1 <- gg.plot(tab, aes(x = hgnc, y = term)) +
    geom_tile(fill = 'gray60', color = 'gray60') +
    geom_point(data = tab.sig, aes(fill = z.show), size = 1.5, pch = 21, show.legend = FALSE) +
    geom_text_repel(data = tab.sig.annot, aes(label = hgnc), size = 3) +
    scale_fill_gradient2(low = 'blue', high = 'red') +
    scale_x_discrete(position = 'top') +
    facet_grid(ontology ~ ., space = 'free', scales = 'free') +
    xlab(cc %&&% ' genes') +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())

.save.png(filename = out.hdr %&&% '_all.png', plot = p1, width = ww, height = hh)
.save.pdf(filename = out.hdr %&&% '_all.pdf', plot = p1, width = ww, height = hh)

p2 <- gg.plot(tab.sig, aes(x = hgnc, y = term)) +
    geom_point(aes(fill = z.show), pch = 22, size = 3) +
    scale_fill_gradient2('medation\neffect\nz-score', low = 'blue', high = 'red') +
    scale_x_discrete(position = 'top') +
    theme(axis.text.x = element_text(angle = 70, hjust = 0, vjust = 0),
          axis.title = element_blank())
          
rr <- tab.sig$category %>% unique() %>% length()
cc <- tab.sig$hgnc %>% unique() %>% length()

ww <- cc * 0.02 + 2 + sapply(tab$term, function(x) nchar(as.character(x))) %>% max() / 10
hh <- rr / 9 + 1

.save.png(filename = out.hdr %&&% '_sig.png', plot = p2, width = ww, height = hh)
.save.pdf(filename = out.hdr %&&% '_sig.pdf', plot = p2, width = ww, height = hh)
