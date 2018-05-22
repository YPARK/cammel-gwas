## confoundeer correction for cis-eqtl calling

LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed 
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)
ldlist := $(shell seq 1 $(NBLOCKS))
TEMPDIR := /broad/hptmp/ypp/cammel-gwas/tempdata

all: $(foreach data, rosmap geuvadis mayo, jobs/step2-$(data)-qtl.txt.gz jobs/step2-$(data)-qtl-mult.txt.gz) jobs/step2-gtex-qtl-mult.txt.gz

long: $(foreach data, rosmap geuvadis mayo, jobs/step2-$(data)-qtl-long.txt.gz jobs/step2-$(data)-qtl-mult-long.txt.gz)

jobs/%-long.txt.gz: jobs/%.txt.gz
	zcat $< | awk 'system(" ! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=36:00:00 -b y -j y -N cis_eqtl_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/step2-geuvadis-qtl.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/geuvadis/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.geuvadis-cis-eqtl.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/geuvadis/" NR "_qtl.txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N geuvadis_cis_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/step2-geuvadis-qtl-mult.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/geuvadis-mult/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.geuvadis-cis-eqtl-mult.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/geuvadis-mult/" NR "_qtl.txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N geuvadis_cis_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/step2-mayo-qtl.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/mayo/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.mayo-cis-eqtl.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/mayo/" NR "_qtl.txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N mayo_cis_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/step2-mayo-qtl-mult.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/mayo-mult/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.mayo-cis-eqtl-mult.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/mayo-mult/" NR "_qtl.txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N mayo_cis_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@


################################################################
jobs/step2-rosmap-qtl.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/rosmap/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.rosmap-cis-eqtl.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/rosmap/" NR "_qtl.txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N rosmap_cis_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/step2-rosmap-qtl-mult.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/rosmap-mult/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.rosmap-cis-eqtl-mult.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/rosmap-mult/" NR "_qtl.txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N rosmap_cis_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
jobs/step2-gtex-qtl-mult.txt.gz: $(LD)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p cis-eqtl/gtex-fqtl-v6/
	@cat $< | tail -n+2 | sed 's/chr//' | awk '{ print "./make.gtex-fqtl-v6.R" FS $$1 FS $$2 FS $$3 FS ("cis-eqtl/gtex-fqtl-v6/" NR "_qtl.txt.gz") }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=1:00:00 -b y -j y -N gtex_eqtl -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

