
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InteractionPoweR <a href='https://dbaranger.github.io/InteractionPoweR'><img src="man/figures/logo.png" align="right" height="139"/></a>

<!-- badges: start -->
<!-- badges: end -->

The `{InteractionPoweR}` package conducts power analyses for regression
models in cross-sectional data sets where the term of interest is an
interaction between two variables, also known as ‘moderation’ analyses.
The package includes functions for simulating data sets, conducting
power analyses, and plotting and interpreting results. Notable package
features include (1) the ability to compute power for interactions
between two continuous variables, (2) effect sizes are all specified as
the cross-sectional Pearson’s correlation, (3) simulations do not assume
that the interacting variables are independent, (4) any variable in the
model, including the outcome, can be binary, and (5) analyses can
incorporate the effects of reliability and skew, both of the interacting
variables, as well as of the outcome variable.

**For more information see [documentation and
examples](https://dbaranger.github.io/InteractionPoweR/articles/articles/InteractionPoweRvignette.html),
and the
[FAQ](https://dbaranger.github.io/InteractionPoweR/articles/articles/CommonQuestions.html).**

We also have a website app which implements the major functions,
available
[**here**](https://dbaranger.github.io/InteractionPoweR/index.html).

Please report bugs, issues, or questions as an [Issue on
Github](https://github.com/dbaranger/InteractionPoweR/issues).

## Installation

You can install InteractionPoweR from github with:

``` r
install.packages("devtools")
devtools::install_github("dbaranger/InteractionPoweR")
```

Sometimes there will be a minor installation error, which can be
resolved by using:

``` r
install.packages("devtools")
devtools::install_github("dbaranger/InteractionPoweR/@HEAD")
```

## Basic Usage

The simplest use-case is when all the input parameters are known. We
know the population-level correlation between our predictors (x1 and x2)
and our outcome, we have a smallest effect size of interest in mind for
our interaction effect size, and our sample size is already set (maybe
we are conducting secondary data analysis). Power can be determined with
a single command. **NB** In all these examples we use 1000 simulations
for speed (`n.iter = 1000`), but for robust results we recommend 10,000
simulations (`n.iter = 10000`).

``` r
library(InteractionPoweR)
library(tictoc)
tic()
test_power<-power_interaction(
  n.iter = 1000,            # number of simulations per unique combination of input parameters
  alpha = 0.05,             # alpha, for the power analysis
  N = 350,                  # sample size
  r.x1x2.y = .15,           # interaction effect to test (correlation between x1*x2 and y)
  r.x1.y = .2,              # correlation between x1 and y
  r.x2.y = .1,              # correlation between x2 and y
  r.x1.x2 = .2,             # correlation between x1 and x2
  seed = 581827             # seed, for reproducibility - this generally should not be set
)
#> [1] "Checking for errors in inputs..."
#> [1] "Performing 1000 simulations"
toc()
#> 16.66 sec elapsed
test_power
#> # A tibble: 1 x 1
#>     pwr
#>   <dbl>
#> 1 0.809
```

We see that we have 80.9% power to detect the effect of interest.
