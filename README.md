Annotation  
=======

This repository contains python and R utilities for analysis and downstream visualization of genetic data: 

* parse_manifest
* parse_vep_output
* worldFreq 

## `parse_manifest`
This is a utility for parsing the MEGA manifest provided by Illumina (https://www.illumina.com/science/consortia/human-consortia/multi-ethnic-genotyping-consortium.html). The MEGA platform enables high imputation accuracy for multi-ethnic cohorts and was developed by the PAGE Study (www.pagestudy.org). To read about trans and multi-ethnic GWAS using the MEGA platform, check out our preprint on biorxiv: (http://www.biorxiv.org/cgi/content/short/188094v1). 

## `parse_vep_output`
Script for selecting annotations of interest from a VEP VCF file and simplifying into a tab-delimited TXT file. 
`Usage: parse_vep_output.py <VCF file> <annotation type>`

## `worldFreq`
Analysis and visualization of minor allele frequency or carrier frequency for a large number of globally representative populations. `worldFreq` was used for studying COL27A1.pG697R carrier frequencies in Figure 4A in Belbin GM, Odgis J, Sorokin EP "Genetic identification of a common collagen disease in Puerto Ricans via identity-by-descent mapping in a health system" (2017). (https://elifesciences.org/articles/25060/figures) And also for studying HLA-B star 57:01 haplotype diversity  Figure 5 in Wojcik et al, "Genetic diversity turns a new PAGE in our understanding of complex traits" (2017). (http://www.biorxiv.org/content/early/2017/09/15/188094) 
