# Contents 

Data and scripts of "Natural variation in copper tolerance in *Drosophila melanogaster* is shaped by transcriptional and physiological changes in the gut" (Green et al. 2021).

**Differential gene expression analyses**
1. `DE_analysis.R`: R script for performing DE analyses, GO enrichment and visualizations
2.  **DESeq2-input** folder: `sample_sheet_*.txt` contain the table of sample information for DESeq2
3.  **DESeq2-input** folder: `counts_*.txt` are the raw gene counts 

**Linear model analyses**
1. `analysisLM.Rmd`: Rmarkdown for performing the multiple linear regression model 
2. `dataLM.tab` data for performing the multiple linear regression model 

**Kaplan-Meier analyses**
1. **Kaplan-Meier** folder: `data_to_kaplan.pl` script to create input for survival analysis
2. **Kaplan-Meier** folder: `survival.R` R script to analyse mortality data (Kaplan-Meier plots and Log Rank Test in R)

# Reference

Green et al., (2021) The genomic basis of copper tolerance in *Drosophila* is shaped by a complex interplay of regulatory and environmental factors ([BioRxiv](https://www.biorxiv.org/content/10.1101/2021.07.12.452058v1)). 
