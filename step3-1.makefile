## Test trans-effect mediated by cis-effect

CUTOFF := 0.01

cisTransFiles := $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, mediation.trans/gene_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).pairs.gz))))

all:

cisTrans: $(cisTransFiles)

Trans: $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, jobs/step3-1_trans_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).txt.gz))))

trans-long: $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, jobs/step3-1_trans_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig)-long.txt.gz))))

Med: $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, jobs/step3-1_mediation_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).txt.gz))))

med-long: $(foreach gwas, igap, $(foreach qtl, rosmap, $(foreach gam, 4, $(foreach eig, 2, jobs/step3-1_mediation_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig)-long.txt.gz))))

################################################################
## 1. Idetify trans pairs
mediation.trans/%.pairs.gz: mediation/%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh Rscript --vanilla ./make.cis-trans.R $< $(CUTOFF) $@

################################################################
## 2. Call trans eQTLs
jobs/step3-1_trans_%-long.txt.gz: jobs/step3-1_trans_%.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=36:00:00 -b y -j y -N trans_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-1_trans_%.txt.gz: mediation.trans/gene_%.pairs.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | sed 's/:/\t/g' | awk '{ print "./make.rosmap-trans-eqtl.R" FS $$1 FS $$2 FS $$3 FS $$4 FS $$5 FS $$6 FS ("trans-eqtl/$*/chr" $$1 "/" $$7 "_" $$8 ".txt.gz") }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N rosmap_trans_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
## 3. Trans mediation

# % = igap_rosmap_gammax-4_eigen-2
jobs/step3-1_mediation_%.txt.gz: mediation.trans/gene_%.pairs.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | sed 's/:/\t/g' | awk -vGWAS=$(shell echo $* | awk -F'_' '{ print $$1 }') -vQTL=$(shell echo $* | awk -F'_' '{ print $$2 }') -vGAM=$(shell echo $* | awk -F'_' '{ print $$3 }' | sed 's/gammax-//g') -vGENO=$(shell echo $* | awk -F'_' '{ print toupper($$2) "_GENO" }') -vEIG=$(shell echo $* | awk -F'_' '{ print $$4 }' | sed 's/eigen-//g') '{ print ("./make.cammel-" GWAS "-trans.R") FS ("trans-eqtl/$*/chr" $$1 "/" $$7 "_" $$8 ".txt.gz") FS GENO FS ("1e" GAM) FS ("1e-" EIG) FS ("mediation.trans/" GWAS "/" QTL "/chr" $$1 "/" $$7 "_" $$8) }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=4:00:00 -b y -j y -N rosmap_trans_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-1_mediation_%-long.txt.gz: jobs/step3-1_mediation_%.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".trans-mediation.txt.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=16:00:00 -b y -j y -N trans_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
## 4. check cis-trans coherence

