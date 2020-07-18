# dbh

An R package for the dependent-adjusted Benjamini-Hochberg and step-up procedures

<!-- badges: start -->
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/lihualei71/dbh?branch=master&svg=true)](https://ci.appveyor.com/project/lihualei71/dbh)
<!-- badges: end -->

## Overview

This R package implements the dependent-adjusted Benjamini-Hochberg procedure, Benjamini-Yekutieli procedure, and step-up procedures for parametric and non-parametric multiple testing problems with exact false discovery rate control. The procedures are proposed in our paper: [Conditional calibration for false discovery rate control under dependence](https://arxiv.org/abs/). The current version supports one- and two-sided tests on  multivariate z-statistics (`dBH_mvgauss`), multivariate t-statistis (`dBH_mvt`), and fixed-design homoscedastic Gaussian linear models (`dBH_lm`). The procedures for other problems will be added in the future versions.

## Installation

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("lihualei71/cfcausal")
```

## Usage Examples

We illustrate the usage of `dBH_mvgauss` below. For details please read the manual.

``` r
library(dbh)

## basic example code
# Generate mu and Sigma for an AR process
n <- 100
rho <- 0.8
Sigma <- rho^(abs(outer(1:n, 1:n, "-")))
mu1 <- 2.5
nalt <- 10
mu <- c(rep(mu1, nalt), rep(0, n - nalt))

# Generate the z-values
set.seed(1)
zvals <- rep(NA, n)
zvals[1] <- rnorm(1)
for (i in 2:n){
    zvals[i] <- zvals[i - 1] * rho + rnorm(1) * sqrt(1 - rho^2)
}
zvals <- zvals + mu

# Run dBH_1(\alpha) for one-sided tests
alpha <- 0.05
res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
                   gamma = 1, niter = 1, avals_type = "BH")

# Run dBH_1(\alpha) for two-sided tests
res <- dBH_mvgauss(zvals = zvals, Sigma = Sigma, side = "right", alpha = alpha,
                   gamma = 1, niter = 1, avals_type = "BH") 
```

