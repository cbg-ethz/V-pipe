BIBS = $(wildcard *.bib)
HTMLS = $(patsubst %.bib, ../_includes/literature-%.html, $(BIBS))

TMP = primary.xml
EXT = bib ris end word.xml
DOWN_BASE = ../assets/cite/v-pipe
DOWNLOADS = $(patsubst %, $(DOWN_BASE).%, $(EXT))

.phony: test all clean

all: $(HTMLS) $(DOWNLOADS)

clean:
	rm -v $(HTMLS) $(TMP) $(DOWNLOADS)

$(HTMLS): ref.md

../_includes/literature-%.html:%.bib
	pandoc --citeproc $(patsubst %.csl,--csl=%.csl,$(^:%.bib=--bibliography %.bib)) -o $@

%.xml:%.bib
	bib2xml $^ > $@

$(DOWN_BASE).%:$(TMP)
	xml2$(patsubst $(DOWN_BASE).%,%,$@) $^ > $@

$(DOWN_BASE).word.xml:$(TMP)
	xml2wordbib $^ > $@
