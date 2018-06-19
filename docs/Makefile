RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)

all: $(HTML)

%.html: %.rmd
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

