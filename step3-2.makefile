LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)

all:

run-jobs: $(foreach gwas, igap, jobs/step3-2-$(gwas).txt.gz)

long: $(foreach gwas, igap, jobs/step3-2-$(gwas).long.gz)

jobs/step3-2-%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-igap-joint.R" FS $$1 FS ("1e"g) FS ("1e-"e) FS ("mediation.joint/$*/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=8:00:00 -b y -j y -N cammel_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-2-%.long.gz: jobs/step3-2-%.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".null.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=36:00:00 -b y -j y -N $*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

