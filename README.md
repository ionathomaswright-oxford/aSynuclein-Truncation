# ***⍺-synuclein truncation is a maturity-dependent phenotype in human dopaminergic neurons***

**This repository contains the scripts and metadata necessary to replicate the transcriptomics analyses included in the manuscript.**


---


## Data / code availability

- Raw (FASTQ) and processed (Salmon quant.sf) bulk RNA-seq data: available from the Gene Expression Omnibus (GEO) under accession number GSE353879.
- Metadata: included in this GitHub repository ('XXX'), within folder 'metadata'.
- Scripts: included in this GitHub repository ('XXX'), within folder 'scripts'.
- Software details: included in this GitHub repository ('XXX'), within folder 'software_details'.


---


## Code execution steps

**RWM dataset**
1. rwm/make_feature_metadata.R
2. rwm/salmon.R

**FOUNDIN-PD dataset**
1. foundinpd/make_sample_metadata.R
2. foundinpd/make_feature_metadata.R
3. foundinpd/salmon.R
4. foundinpd/dgea.R


---


## Requirements

- R package `rnaseqtools` available from [GitHub](https://github.com/benvallin/rnaseqtools) with:

```r
remotes::install_github("benvallin/rnaseqtools")
```


---