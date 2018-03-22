## run simulations

LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)

jobs/step4-simulation.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(LD) | tail -n+2  | awk '($$3 - $$2) > 1e6 { print NR }' | awk 'BEGIN{ srand(9689) } { print rand() FS $$1 }' | sort -k1g | head -n20 | awk '{ print $$2 }' | awk '{ hq = 0.17; split("10 50 100", ngenes); split("1 5 10", ncgenes); for(g = 1; g <= 2; ++g) { ng = ngenes[g]; nc = ncgenes[g]; for(sr = 0.3; sr <= 0.6; sr += 0.3) for(gd = 2; gd <= 3; gd ++) for(gg = 1; gg <= 10; gg+=3) for(qq = 2; qq < 5; qq+=2) for(hm = 0.0; hm <= 0.3; hm += 0.05) for(hu = 0.0; hu <= 0.5; hu += 0.1) for(nu = 1; nu <= 3; ++nu) print "./make.simulation.R " $$1 FS sr FS ng FS nc FS qq FS nu FS hq FS hm FS hu FS ("./simulation/" ng "_" nc "/" (ng "_" nc "_" sr "_" nu "_" qq "_" hq "_" hm "_" hu) "_" $$1 ".txt.gz") } }' | awk 'system("! [ -f " $$NF " ]") == 0' | head # gzip > $@
##	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_simulation -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step4-simulation-long.txt.gz: jobs/step4-simulation.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=36:00:00 -b y -j y -N cammel_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@
