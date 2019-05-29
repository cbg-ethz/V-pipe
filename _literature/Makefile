BIBS = $(wildcard *.bib)
HTMLS = $(patsubst %.bib, ../_includes/literature-%.html, $(BIBS))

.phony: test all

all: $(HTMLS)

$(HTMLS): ref.md

../_includes/literature-%.html:%.bib
	pandoc --filter=pandoc-citeproc $(patsubst %.csl,--csl=%.csl,$(^:%.bib=--bibliography %.bib)) -o $@

