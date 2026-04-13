# TPSAS: Two-Piece Sinh-Arcsinh Distributions (R package)

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`TPSAS` is an R package implementing the family of **two-piece sinh-arcsinh
(TPSAS) distributions**, introduced in:

> Rubio, F.J., Ogundimu, E.O., and Hutton, J.L. (2016). On modelling asymmetric
> data using two-piece sinh–arcsinh distributions. *Brazilian Journal of
> Probability and Statistics* **30**(4): 485–501.
> [DOI: 10.1214/15-BJPS290](https://doi.org/10.1214/15-BJPS290)

The TPSAS family combines the **two-piece construction**, which introduces
asymmetry by joining two rescaled half-densities at the mode, with the (symmetric)
**sinh-arcsinh transformation** of Jones and Pewsey (2009), which controls the tail weight through a shape parameter. The
result is a flexible four-parameter family of unimodal distributions on
$\mathbb{R}$ capable of accommodating a wide range of asymmetry and tail
behaviour, making it well suited for modelling skewed continuous data.

## Installation

```r
# install.packages("devtools")
devtools::install_github("FJRubio67/TPSAS")
library(TPSAS)
```

## Main functions

The package provides density (`d`), distribution (`p`), quantile (`q`), and
random sampling (`r`) functions for the TPSAS family:

| Function | Description |
|---|---|
| `dtpsas` | Probability density function |
| `ptpsas` | Cumulative distribution function |
| `qtpsas` | Quantile function |
| `rtpsas` | Random number generation |

The distribution is parameterised by location $\mu$, scale $\sigma > 0$,
skewness $\varepsilon$, and tail weight $\delta > 0$. Setting
$\varepsilon = 0$ and $\delta = 1$ recovers the symmetric sinh-arcsinh
baseline. For full documentation: `?dtpsas`

## Quick example

```r
library(TPSAS)

# Density with positive skewness and heavy tails
x <- seq(-5, 8, length.out = 500)
plot(x, dtpsas(x, mu = 0, sigma = 1, epsilon = 0.5, delta = 0.8),
     type = "l", ylab = "Density", main = "TPSAS distribution")

# Random sampling
sim <- rtpsas(1000, mu = 0, sigma = 1, epsilon = 0.5, delta = 0.8)
```

## Tutorials

- [TPSAS package](https://rpubs.com/FJRubio/TPSAS) — illustrative walkthrough
  with examples and parameter interpretation
- [The original sinh-arcsinh distribution](https://rpubs.com/FJRubio/SAS) —
  background on the Jones-Pewsey baseline distribution

## Citation

If you use this package, please cite:

```bibtex
@article{rubio:2016,
  author  = {Rubio, F.J. and Ogundimu, E.O. and Hutton, J.L.},
  title   = {On modelling asymmetric data using two-piece sinh--arcsinh
             distributions},
  journal = {Brazilian Journal of Probability and Statistics},
  volume  = {30},
  number  = {4},
  pages   = {485--501},
  year    = {2016},
  doi     = {10.1214/15-BJPS290}
}
```

## Related packages

- [twopiece](https://github.com/FJRubio67/twopiece) — the broader family of
  two-piece distributions (parent package)
- [DTP](https://github.com/FJRubio67/DTP) — double two-piece distributions,
  including the double two-piece sinh-arcsinh
- [FlexBinReg](https://github.com/FJRubio67/FlexBinReg) — binary regression
  using two-piece link functions

## License

This package is licensed under the [MIT License](LICENSE).
