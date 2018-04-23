## Test trans-effect mediated by cis-effect

CUTOFF := 0.01

cisTransFiles := $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, mediation.trans/gene_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).pairs.gz))))

all:

cisTrans: $(cisTransFiles)

Trans: $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, jobs/step3-1_trans_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).txt.gz))))

################################################################
## 1. Idetify trans pairs
mediation.trans/%.pairs.gz: mediation/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh Rscript --vanilla ./make.cis-trans.R $< $(CUTOFF) $@

################################################################
## 2. Call trans eQTLs

jobs/step3-1_trans_%.txt.gz: mediation.trans/gene_%.pairs.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | sed 's/:/\t/g' | awk '{ print "./make.rosmap-trans-eqtl.R" FS $$1 FS $$2 FS $$3 FS $$4 FS $$5 FS $$6 FS ("trans-eqtl/$*/chr" $$1 "/" $$7 "_" $$8 ".txt.gz") }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N rosmap_trans_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@
