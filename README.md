GENESIS: GENetic Effect-Size distribution Inference from Summary-level data
====
  
## Overview

The goal of GENESIS is to analyze summary-level GWAS statistics and external linkage disequilibrium information to estimate common variants effect-size distributions, characterized by the proportion of underlying susceptibility SNPs and a flexible normal-mixture model for their effects. This package allows flexibility by considering a 2- or 3-component model, which respectively incorporate a single normal, or a mixture of two normal distributions, for specifying the effects of non-null SNPs. This package also allows users to make predictions regarding yield of future GWAS with larger sample sizes.

## GENESIS Installation

GENESIS software can be installed via Github. To install the latest version of GENESIS package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("yandorazhang/GENESIS")
```

## Quality Control
The input GWAS summary statistics should contain 3 columns: 

1. SNP rsID, 
2. original z-statistics got from GWAS study (uncorrected by genomic control factor), 
3. effective sample size of GWAS study, which can vary for different SNPs; for disease traits, the sample size should be effective sample size, i.e., (# of cases)*(# of controls)/(total # of cases & controls).

The input GWAS summary statistics are strongly recommended to do filtering before fitting to the model: 

1. If sample size is different for different SNPs, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size. 
2. Filter SNPs to Hapmap3 SNPs which are not in the major histocompatibility complex (MHC) region. For Hapmap3 SNP list without MHC region, type ```data(w_hm3.noMHC.snplist)``` in ```R```.
3. Remove SNPs with extremely large effect sizes (z^2 > 80).




## Citation

Please cite the following paper when you use GENESIS:


[Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits." Nature genetics 50.9 (2018): 1318.] (https://www.nature.com/articles/s41588-018-0193-x)



## Contact the Author
Software Developer/Maintainer: Yan Zhang (yzhan284@jhu.edu or  yandorazhang@gmail.com)


