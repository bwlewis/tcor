# tcor
An R package for fast and memory-efficient computation of thresholded correlation matrices

A preprint of the companion paper is available here: http://arxiv.org/abs/1512.07246


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
