RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)

all: $(HTML) $(PDF)

%.html: %.rmd
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

