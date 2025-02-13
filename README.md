# recall (Calibrated Clustering with Artificial Variables) Reproducibility <img src="man/figures/recall_logo.png" align="right" alt="" width="120"/>


## Introduction

We hope to make it as simple as possible to reproduce the results found in our paper `Artificial variables help to avoid over-clustering in single-cell RNA-sequencing`.
To that end, we have organized our analysis scripts as a series of `R` vignettes in this repository.

## Reproducing Figures from the `recall` manuscript

Clone the repository:

```bash
git clone https://github.com/lcrawlab/recallreproducibility
```

Then, navigate to the repo directory and launch R
```bash
cd recallreproducibility
```

```bash
R
```

You can build the entire website using the following line of `R` code.
```r
pkgdown::build_site()
```

Note that each `Rmarkdown` file is not fully run by `R` by default. To properly run an an `Rmarkdown` file run during the website building process, you need to remove 

```R
knitr::opts_chunk$set(eval = FALSE)
```

from the file header. We do not recommend doing this for the vignettes that actually use `recall`, `sc-SHC`, and `CHOIR` for clustering because they have a long runtime. Rather, the R portions of these files should be put in a script and run using `Rscript`.


## Relevant Citations
`recall` is currently on the bioRxiv, [here](https://www.biorxiv.org/content/10.1101/2024.03.08.584180v1).

A. DenAdel, M. Ramseier, A. Navia, A. Shalek, S. Raghavan, P. Winter, A. Amini, and L. Crawford. A knockoff calibration method to avoid over-clustering in single-cell RNA-sequencing. _bioRxiv_.

## Questions and Feedback
For questions or concerns with `recallreproducibility` or the `recall` `R` package, please contact
[Alan DenAdel](mailto:alan_denadel@brown.edu) or [Lorin Crawford](lcrawford@microsoft.com). Any feedback on the manuscript or figure reproducibility is greatly appreciated.
