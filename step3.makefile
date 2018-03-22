## run CaMMEL

LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)

all:

prepare:
	qsub -P compbio_lab -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N ad_gwas ./run.sh ./make.distribute-ad.R

run-jobs: $(foreach data, mayo rosmap geuvadis gtex, jobs/step3-ad-$(data).txt.gz)

run-long: $(foreach data, mayo rosmap geuvadis gtex, jobs/step3-ad-$(data)-long.txt.gz)

figure-jobs: $(foreach data, rosmap mayo geuvadis, jobs/step3-ad-$(data)_show.txt.gz)

table: $(foreach qtl, mayo rosmap geuvadis GTExFQTLv6, $(foreach gam, 4, $(foreach eig, 2, mediation/gene_ad_$(qtl)_gammax-$(gam)_eigen-$(eig).txt.gz)))

################################################################
jobs/%-long.txt.gz: jobs/%.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".mediation.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=36:00:00 -b y -j y -N cammel_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-rosmap.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/rosmap/" $$1 "_qtl.txt.gz" FS 356 FS "ROSMAP_GENO" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/rosmap/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_rosmap_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-geuvadis.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/geuvadis/" $$1 "_qtl.txt.gz" FS 358 FS "1KG_EUR" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/geuvadis/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_geuvadis_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-mayo.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/mayo/" $$1 "_qtl.txt.gz" FS 266 FS "Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/mayo/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_mayo_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-rosmap-mult.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/rosmap-mult/" $$1 "_qtl.txt.gz" FS 356 FS "ROSMAP_GENO" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/rosmapMult/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_rosmap_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-geuvadis-mult.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/geuvadis-mult/" $$1 "_qtl.txt.gz" FS 358 FS "1KG_EUR" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/geuvadisMult/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_geuvadis_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-mayo-mult.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/mayo-mult/" $$1 "_qtl.txt.gz" FS 266 FS "Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/mayoMult/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_mayo_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad-gtex.txt.gz:
	echo "nothing to do" | gzip > $@

jobs/step3-ad-gtex-mult.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-ad.R" FS $$1 FS "cis-eqtl/gtex-fqtl-v6/" $$1 "_qtl.txt.gz" FS 450 FS "GTEX_GENO" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/GTExFQTLv6/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_gtex_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

################################################################
mediation/gene_%.txt.gz: make.cammel-gene-table.R
	./run.sh Rscript --vanilla $< mediation/$(shell echo $* | sed 's/_/\//g' | sed 's/-/_/g' | sed 's/\./-/g') $@

################################################################
# visualize significant LD blocks
FSR := 0.5

jobs/step3-ad-%_show.txt.gz: $(foreach g, 4, $(foreach e, 2, jobs/step3-ad_%.gammax-$(g)_eigen-$(e).show.txt.gz))
	@cat $^ > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=8:00:00 -b y -j y -N figure_ad_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-ad_rosmap.%.show.txt.gz: mediation/gene_ad_rosmap_%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@tail -n+2 $(LD) | awk '{ print "_ld" "\t" NR "\t" $$1 "\t" $$2 "\t" $$3 }' | gzip | zcat - $< | awk -v args="$(shell echo $* | awk -F'_' '{ gsub("gammax-","", $$1); gsub("eigen-","", $$2); print ("1e" $$1) "\t" ("1e-" $$2) }')" -F'\t' '($$1 == "_ld") { k = ($$3 FS $$4 FS $$5); ld[k] = $$2; } ($$1 != "_ld") && /igap/ && ($$(NF -1) < $(FSR)) { _k = "chr" $$25 FS $$26 FS $$27; k = ld[_k]; show[k]++; } END { for(k in show) print "./show.cammel-ad-local.R" FS k FS ("cis-eqtl/rosmap/" k "_qtl.txt.gz") FS 356 FS "ROSMAP_GENO" FS args FS "$<" FS ("figure/gene/ad/rosmap/$*/Fig_Loc")  }' | sort -k2n | awk 'system("[ -f " $$(NF - 1) " ]") == 0' | gzip > $@

jobs/step3-ad_mayo.%.show.txt.gz: mediation/gene_ad_mayo_%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@tail -n+2 $(LD) | awk '{ print "_ld" "\t" NR "\t" $$1 "\t" $$2 "\t" $$3 }' | gzip | zcat - $< | awk -v args="$(shell echo $* | awk -F'_' '{ gsub("gammax-","", $$1); gsub("eigen-","", $$2); print ("1e" $$1) "\t" ("1e-" $$2) }')" -F'\t' '($$1 == "_ld") { k = ($$3 FS $$4 FS $$5); ld[k] = $$2; } ($$1 != "_ld") && /igap/ && ($$(NF -1) < $(FSR)) { _k = "chr" $$25 FS $$26 FS $$27; k = ld[_k]; show[k]++; } END { for(k in show) print "./show.cammel-ad-local.R" FS k FS ("cis-eqtl/mayo/" k "_qtl.txt.gz") FS 266 FS "Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed" FS args FS "$<" FS ("figure/gene/ad/mayo/$*/Fig_Loc")  }' | sort -k2n | awk 'system("[ -f " $$(NF - 1) " ]") == 0' | gzip > $@

jobs/step3-ad_geuvadis.%.show.txt.gz: mediation/gene_ad_geuvadis_%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@tail -n+2 $(LD) | awk '{ print "_ld" "\t" NR "\t" $$1 "\t" $$2 "\t" $$3 }' | gzip | zcat - $< | awk -v args="$(shell echo $* | awk -F'_' '{ gsub("gammax-","", $$1); gsub("eigen-","", $$2); print ("1e" $$1) "\t" ("1e-" $$2) }')" -F'\t' '($$1 == "_ld") { k = ($$3 FS $$4 FS $$5); ld[k] = $$2; } ($$1 != "_ld") && /igap/ && ($$(NF -1) < $(FSR)) { _k = "chr" $$25 FS $$26 FS $$27; k = ld[_k]; show[k]++; } END { for(k in show) print "./show.cammel-ad-local.R" FS k FS ("cis-eqtl/geuvadis/" k "_qtl.txt.gz") FS 358 FS "1KG_EUR" FS args FS "$<" FS ("figure/gene/ad/geuvadis/$*/Fig_Loc")  }' | sort -k2n | awk 'system("[ -f " $$(NF - 1) " ]") == 0' | gzip > $@

