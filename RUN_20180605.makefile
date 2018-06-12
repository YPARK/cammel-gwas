LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)

all:

prepare:
	qsub -P compbio_lab -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N igap_gwas ./run.sh ./make.distribute-gwas.R
	qsub -P compbio_lab -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N igap_gwas ./run.sh ./make.distribute-gwas-ptsd.R

GWAS := pgc ukbb igap

run-jobs: $(foreach gwas, $(GWAS), jobs/20180605-$(gwas).txt.gz)

run-null: $(foreach gwas, $(GWAS), jobs/20180605-$(gwas).null.gz)

long: $(foreach gwas, $(GWAS), jobs/20180605-$(gwas).long.gz)

jobs/20180605-%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-$*-joint.R" FS $$1 FS ("1e"g) FS ("1e-"e) FS "rosmap-mult:mayo-mult:geuvadis-mult" FS ("result/20180605/obs/$*/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180605-%.null.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-$*-joint.R" FS $$1 FS ("1e"g) FS ("1e-"e) FS "rosmap-mult:mayo-mult:geuvadis-mult" FS ("result/20180605/null/$*/" $$1) FS "TRUE" }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N null_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180605-igap.long.gz: jobs/20180605-igap.txt.gz
	zcat $< | awk 'system("! [ -f " $$6 ".igap_ad.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=72:00:00 -b y -j y -N igap_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180605-pgc.long.gz: jobs/20180605-pgc.txt.gz
	zcat $< | awk 'system("! [ -f " $$6 ".pgc_ptsd_all_ea.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=72:00:00 -b y -j y -N pgc_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/20180605-ukbb.long.gz: jobs/20180605-ukbb.txt.gz
	zcat $< | awk 'system("! [ -f " $$6 ".ukbb_moodswings.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=72:00:00 -b y -j y -N ukbb_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
table: $(foreach _trait, ukbb_ad_mat ukbb_ad_pat ukbb_moodswings ukbb_neuroticism pgc_scz pgc_mdd pgc_asd pgc_adhd pgc_ocd pgc_ptsd_civ_ea igap_ad, result/20180605/obs/gene-$(_trait).txt.gz)

# % = igap-gammax_4-eigen_2.mediation

result/20180605/obs/gene-%.txt.gz: make.cammel-gene-table.R
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh Rscript --vanilla $< result/20180605/obs/$(shell echo $* | awk -F'_' '{ print $$1 "/ " $$0 }') $@

