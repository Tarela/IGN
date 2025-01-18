# IGN: invariable gene set based normalization for chromatin accessibility profile data analysis

Chromatin accessibility profiles generated using ATAC-seq or DNase-seq carry functional information of the regulatory genome that control gene expression. Appropriate normalization of ATACseq and DNase-seq data is essential for accurate differential analysis when studying chromatin dynamics. Existing normalization methods usually assume the same distribution of genomic signals across samples; however, this approach may not be appropriate when there are global changes of chromatin accessibility levels between experimental conditions/samples.

We present IGN (Invariable Gene Normalization), a method for ATAC-seq and DNase-seq data normalization. IGN is performed by normalizing the promoter chromatin accessibility signals for a set of genes that are unchanged in expression, usually obtained from accompanying RNA-seq data, and extrapolating to normalize the genome-wide chromatin accessibility profile. We demonstrate that IGN outperforms existing methods. As the first chromatin accessibility normalization method allowing for global difference, IGN can be widely applied for differential ATAC-seq and DNase-seq analysis.

[![R 3.0](https://img.shields.io/badge/R-3.0-blue.svg)](https://www.r-project.org/)
[![GitHub license](https://img.shields.io/github/license/Tarela/IGN)](https://github.com/Tarela/IGN/blob/main/LICENSE)

## 0. Introduction of IGN package
IGN is performed by normalizing the promoter chromatin accessibility signals for a given gene set that is unchanged in expression, usually obtained from accompanying RNA-seq data, and extrapolating to scale the genome-wide chromatin accessibility profile. This function allows users to normalize ATAC/DNase-seq signal matrix based on promoter signal of invariable genes. 

- Changelog<br>
v1.0.0 First version of IGN.
v1.1.0 CSBJ revision
## 1. Installation
- Package requirements<br>
IGN requires [Rscript](https://www.r-project.org) v3+ to run.<br>

- Github installation
Simply clone this repo and install the IGN package in R
```sh
$ git clone https://github.com/Tarela/IGN.git
$ cd IGN/
$ R CMD INSTALL IGN_1.0.0.tar.gz
```

## 2. Run IGN (usage)
- Check the R help document in the IGN package for detailed tutorial<br>
Type the following code in an R session after installing it:
```R
> library(IGN)
> ?IGN
```
