LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)

all:

prepare:
	qsub -P compbio_lab -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N igap_gwas ./run.sh ./make.distribute-gwas.R
	qsub -P compbio_lab -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N igap_gwas ./run.sh ./make.distribute-gwas-ptsd.R

run-jobs: $(foreach gwas, pgc igap ukbb, jobs/step3-$(gwas).txt.gz)

long: $(foreach gwas, pgc igap ukbb, jobs/step3-$(gwas).long.gz)

jobs/step3-%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; for(e = 3; e >= 2; --e) print "./make.cammel-$*-joint.R" FS $$1 FS ("1e"g) FS ("1e-"e) FS ("mediation.joint/$*/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=8:00:00 -b y -j y -N cammel_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-igap.long.gz: jobs/step3-igap.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".null.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=72:00:00 -b y -j y -N igap_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-igap.long.gz: jobs/step3-pgc.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".pgc_mdd.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=72:00:00 -b y -j y -N igap_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ukbb.long.gz: jobs/step3-ukbb.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".moodswings.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=72:00:00 -b y -j y -N ukbb_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

