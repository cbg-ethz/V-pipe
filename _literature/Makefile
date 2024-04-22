BIBS = $(wildcard *.bib)
HTMLS = $(patsubst %.bib, ../_includes/literature-%.html, $(BIBS))

.phony: test all clean

all: $(HTMLS)

clean:
	rm -v $(HTMLS)

$(HTMLS): ref.md

../_includes/literature-%.html:%.bib
	pandoc --citeproc $(patsubst %.csl,--csl=%.csl,$(^:%.bib=--bibliography %.bib)) -o $@

