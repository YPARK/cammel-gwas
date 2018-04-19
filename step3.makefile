## run CaMMEL

LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)

all:

prepare:
	qsub -P compbio_lab -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N igap_gwas ./run.sh ./make.distribute-ad.R

run-jobs: $(foreach gwas, ptsdEA igap, $(foreach data, mayo rosmap, jobs/step3-$(gwas)-$(data).txt.gz))

run-long: $(foreach gwas, ptsdEA igap, $(foreach data, mayo rosmap, jobs/step3-$(gwas)-$(data)-long.txt.gz))

figure-jobs: $(foreach data, rosmap mayo, jobs/step3-igap-$(data)_show.txt.gz)

table: $(foreach gwas, ptsdEA igap, $(foreach qtl, mayo rosmap, $(foreach gam, 4, $(foreach eig, 2, mediation/gene_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig).txt.gz))))

GO: $(foreach gwas, igap, $(foreach qtl, mayo rosmap, $(foreach gam, 4, $(foreach eig, 2, mediation/gene_$(gwas)_$(qtl)_gammax-$(gam)_eigen-$(eig)_go_tab.txt.gz))))

mediation/%_go_tab.txt.gz: mediation/%.txt.gz
	./run.sh Rscript --vanilla $< $(shell echo $@ | sed 's/_tab.txt.gz//')

################################################################
jobs/%-long.txt.gz: jobs/%.txt.gz
	zcat $< | awk 'system("! [ -f " $$NF ".mediation.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=36:00:00 -b y -j y -N cammel_$*_long -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-%-rosmap.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-$*.R" FS $$1 FS "cis-eqtl/rosmap/" $$1 "_qtl.txt.gz" FS "ROSMAP_GENO" FS ("1e"g) FS ("1e-"e) FS ("mediation/$*/rosmap/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_rosmap_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-%-geuvadis.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-$*.R" FS $$1 FS "cis-eqtl/geuvadis/" $$1 "_qtl.txt.gz" FS "1KG_EUR" FS ("1e"g) FS ("1e-"e) FS ("mediation/$*/geuvadis/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_geuvadis_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-%-mayo.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-$*.R" FS $$1 FS "cis-eqtl/mayo/" $$1 "_qtl.txt.gz" FS "Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed" FS ("1e"g) FS ("1e-"e) FS ("mediation/$*/mayo/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_mayo_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-%-gtex.txt.gz:
	echo "nothing to do" | gzip > $@

################################################################
mediation/gene_%.txt.gz: make.cammel-gene-table.R
	./run.sh Rscript --vanilla $< mediation/$(shell echo $* | sed 's/_/\//g' | sed 's/-/_/g' | sed 's/\./-/g') $@

################################################################
# visualize significant LD blocks
FSR := 0.1

jobs/step3-igap-%_show.txt.gz: $(foreach g, 4, $(foreach e, 2, jobs/step3-igap_%.gammax-$(g)_eigen-$(e).show.txt.gz))
	@cat $^ > $@
	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=8:00:00 -b y -j y -N figure_igap_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

jobs/step3-igap_rosmap.%.show.txt.gz: mediation/gene_igap_rosmap_%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@tail -n+2 $(LD) | awk '{ print "_ld" "\t" NR "\t" $$1 "\t" $$2 "\t" $$3 }' | gzip | zcat - $< | awk -v args="$(shell echo $* | awk -F'_' '{ gsub("gammax-","", $$1); gsub("eigen-","", $$2); print ("1e" $$1) "\t" ("1e-" $$2) }')" -F'\t' '($$1 == "_ld") { k = ($$3 FS $$4 FS $$5); ld[k] = $$2; } ($$1 != "_ld") && /IGAP/ && ($$(NF -1) < $(FSR)) { _k = "chr" $$24 FS $$25 FS $$26; k = ld[_k]; show[k]++; } END { for(k in show) print "./show.cammel-igap-local.R" FS k FS ("cis-eqtl/rosmap/" k "_qtl.txt.gz") FS 356 FS "ROSMAP_GENO" FS args FS "$<" FS ("figure/gene/ad/rosmap/$*/Fig_Loc")  }' | sort -k2n | gzip > $@

jobs/step3-igap_mayo.%.show.txt.gz: mediation/gene_igap_mayo_%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@tail -n+2 $(LD) | awk '{ print "_ld" "\t" NR "\t" $$1 "\t" $$2 "\t" $$3 }' | gzip | zcat - $< | awk -v args="$(shell echo $* | awk -F'_' '{ gsub("gammax-","", $$1); gsub("eigen-","", $$2); print ("1e" $$1) "\t" ("1e-" $$2) }')" -F'\t' '($$1 == "_ld") { k = ($$3 FS $$4 FS $$5); ld[k] = $$2; } ($$1 != "_ld") && /IGAP/ && ($$(NF -1) < $(FSR)) { _k = "chr" $$24 FS $$25 FS $$26; k = ld[_k]; show[k]++; } END { for(k in show) print "./show.cammel-igap-local.R" FS k FS ("cis-eqtl/mayo/" k "_qtl.txt.gz") FS 266 FS "Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed" FS args FS "$<" FS ("figure/gene/ad/mayo/$*/Fig_Loc")  }' | sort -k2n | gzip > $@


################################################################
# old code
# jobs/step3-igap-rosmap-mult.txt.gz:
# 	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
# 	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-igap.R" FS $$1 FS "cis-eqtl/rosmap-mult/" $$1 "_qtl.txt.gz" FS 356 FS "ROSMAP_GENO" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/rosmapMult/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
# 	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_rosmap_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

# jobs/step3-igap-geuvadis-mult.txt.gz:
# 	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
# 	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-igap.R" FS $$1 FS "cis-eqtl/geuvadis-mult/" $$1 "_qtl.txt.gz" FS 358 FS "1KG_EUR" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/geuvadisMult/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
# 	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_geuvadis_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

# jobs/step3-igap-mayo-mult.txt.gz:
# 	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
# 	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-igap.R" FS $$1 FS "cis-eqtl/mayo-mult/" $$1 "_qtl.txt.gz" FS 266 FS "Mayo/MayoRNAseq_RNAseq_Genome-Wide_Genotypes_HRCimputed" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/mayoMult/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
# 	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_mayo_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

# jobs/step3-igap-gtex-mult.txt.gz:
# 	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
# 	seq 1 $(NBLOCKS) | awk '{ g = 4; e = 2; print "./make.cammel-igap.R" FS $$1 FS "cis-eqtl/gtex-fqtl-v6/" $$1 "_qtl.txt.gz" FS 450 FS "GTEX_GENO" FS ("1e"g) FS ("1e-"e) FS ("mediation/ad/GTExFQTLv6/gammax_" g "/eigen_" e "/" $$1) }' | gzip >> $@
# 	qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=4:00:00 -b y -j y -N cammel_gtex_ad -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

