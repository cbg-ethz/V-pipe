

# Fetching BIB from DOI

https://doi2bib.org/

# Fetching BIB from PubMed

```bash
curl -L "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${PMID}&retmode=xml" | xsltproc pubmed2bibtex.xsl -
``` 

# Convert RIS to BIB

https://www.bruot.org/ris2bib/

# Finally

Do not forget to compile the HTML includes:

```bash
make
```
