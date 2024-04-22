<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version = "1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!--
  - This Extensible Stylesheet Language Transformation file translates XML files
  - as provided by PubMed into BibTeX files.
  -
  - This file was written by Thomas Fischer <fischer@unix-ag.uni-kl.de>
  - It is released under the GNU Public License version 2 or later.
  -
  - To run test this transformation file, run e.g.
  - wget 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=21351276&retmode=xml' -O - | xsltproc  pubmed2bibtex.xsl -
  -->

<xsl:output method="text" omit-xml-declaration="yes" indent="no" encoding="UTF-8"/>
<xsl:strip-space elements="*"/>

<xsl:variable name="smallcase" select="'abcdefghijklmnopqrstuvwxyz'" />
<xsl:variable name="uppercase" select="'ABCDEFGHIJKLMNOPQRSTUVWXYZ'" />


<!-- START HERE -->
<xsl:template match="/">
<!-- process each entry -->
<xsl:apply-templates select="PubmedArticleSet/PubmedArticle" />
</xsl:template>

<xsl:template match="PubmedArticle">
<!-- assuming that there are only journal references -->
<xsl:text>@article{pmid</xsl:text>
<xsl:value-of select="MedlineCitation/PMID" />
<xsl:apply-templates select="MedlineCitation/Article" />
<xsl:apply-templates select="PubmedData/ArticleIdList/ArticleId" />
<xsl:if test="string-length(MedlineCitation/MedlineJournalInfo/NlmUniqueID) > 0"><xsl:text>,
    nlmuniqueid = {</xsl:text><xsl:value-of select="MedlineCitation/MedlineJournalInfo/NlmUniqueID" /><xsl:text>}</xsl:text></xsl:if>
<xsl:text>
}

</xsl:text>
</xsl:template>


<xsl:template match="ArticleId">
<xsl:choose>
<xsl:when test="@IdType='doi'">
<xsl:text>,
    doi = {</xsl:text><xsl:value-of select="." /><xsl:text>}</xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:text>,
    </xsl:text><xsl:value-of select="@IdType" /><xsl:text> = {</xsl:text><xsl:value-of select="." /><xsl:text>}</xsl:text>
</xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="Article">
<xsl:text>,
    title = {{</xsl:text><xsl:value-of select="ArticleTitle" /><xsl:text>}}</xsl:text>
<xsl:apply-templates select="AuthorList" />
<xsl:apply-templates select="Journal" />
<xsl:text>,
    pages = {</xsl:text><xsl:value-of select="Pagination/MedlinePgn" /><xsl:text>}</xsl:text>
<xsl:if test="string-length(Abstract/AbstractText) > 0"><xsl:text>,
    abstract = {</xsl:text><xsl:value-of select="Abstract/AbstractText" /><xsl:text>}</xsl:text></xsl:if>
</xsl:template>



<xsl:template match="Journal">
<!-- going for the journal title's abbreviation, looks better -->
<xsl:text>,
    journal = {{</xsl:text><xsl:value-of select="ISOAbbreviation" /><xsl:text>}}</xsl:text>
<xsl:if test="string-length(JournalIssue/ISSN) > 0"><xsl:text>,
    issn = {</xsl:text><xsl:value-of select="ISSN" /><xsl:text>}</xsl:text></xsl:if>
<xsl:if test="string-length(JournalIssue/Volume) > 0"><xsl:text>,
    volume = {</xsl:text><xsl:value-of select="JournalIssue/Volume" /><xsl:text>}</xsl:text></xsl:if>
<xsl:if test="string-length(JournalIssue/Issue) > 0"><xsl:text>,
    number = {</xsl:text><xsl:value-of select="JournalIssue/Issue" /><xsl:text>}</xsl:text></xsl:if>
<xsl:if test="string-length(JournalIssue/PubDate/Year) > 0"><xsl:text>,
    year = {</xsl:text><xsl:value-of select="JournalIssue/PubDate/Year" /><xsl:text>}</xsl:text></xsl:if>
<xsl:if test="string-length(JournalIssue/PubDate/Month) > 0"><xsl:text>,
    month = </xsl:text><xsl:value-of select="translate(JournalIssue/PubDate/Month, $uppercase, $smallcase)" /><xsl:text></xsl:text></xsl:if>
</xsl:template>


<xsl:template match="AuthorList">
<xsl:text>,
    author = {</xsl:text>
<xsl:apply-templates select="Author"/>
<xsl:if test="@CompleteYN = 'N'"><xsl:text> and others</xsl:text></xsl:if>
<xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="Author">
<xsl:if test="position() > 1"><xsl:text> and </xsl:text></xsl:if>
<xsl:apply-templates select="ForeName"/><xsl:text> </xsl:text>
<xsl:apply-templates select="LastName"/>
</xsl:template>

</xsl:stylesheet>
