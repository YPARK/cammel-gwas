## Test trans-effect mediated by cis-effect

CUTOFF := 0.01

cisTransFiles := $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, mediation.trans/gene_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).pairs.gz))))

all:

cisTrans: $(cisTransFiles)

mediation.trans/%.pairs.gz: mediation/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh Rscript --vanilla ./make.cis-trans.R $< $(CUTOFF) $@



