# tcor
An R package for fast and memory-efficient computation of thresholded correlation matrices

A preprint of the companion note is available from http://arxiv.org/abs/1512.07246
(the article is being revised, we found some cool related references that we're adding...)


## Installation

The package depends on the 
`irlba` (https://cran.r-project.org/web/packages/irlba/)
and  `foreach` (https://cran.r-project.org/web/packages/foreach/)
packages, each available on CRAN.
You can install `tcor` using the `devtools` package
(https://cran.r-project.org/web/packages/devtools/) with:
```r
devtools::install_github("bwlewis/tcor")
```

The algorithm can optionally make use of `foreach` parallel "back-ends." Many
are available, including `doMC`, `doParallel`, and `doRedis`. See the CRAN high
performance computing task view for more info
https://cran.r-project.org/web/views/HighPerformanceComputing.html.

## Example

See the vignette https://github.com/bwlewis/tcor/blob/master/vignettes/brca.Rmd
for an example that uses tcor to compute the most correlated gene expression
vectors from TCGA RNASeq data.
