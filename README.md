
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
model, including the outcome, can have anywhere from 2 (i.e., binary) to
20 discrete values, and (5) analyses can incorporate the effects of
reliability and skew, both of the interacting variables, as well as of
the outcome variable.

**For more information see [documentation and
examples](https://dbaranger.github.io/InteractionPoweR/articles/articles/InteractionPoweRvignette.html),
and the
[FAQ](https://dbaranger.github.io/InteractionPoweR/articles/articles/CommonQuestions.html).**

We also have a website app which implements the major functions,
available [**here**](https://intmoddev.shinyapps.io/intPower/).

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

If you get an error about a corrupt .rdb file, try restarting your R
session.

## Basic Usage

The simplest use-case is when all the input parameters are known. We
know the population-level correlation between our predictors (x1 and x2)
and our outcome, we have a smallest effect size of interest in mind for
our interaction effect size, and our sample size is already set (maybe
we are conducting secondary data analysis). Power can be determined with
a single command.

First - analytic power, using variable correlations (and reliability, if
provided) to estimate how much additional variance is explained by the
interaction term.

``` r
library(InteractionPoweR)

test_power<-power_interaction_r2(
  alpha = 0.05,             # alpha, for the power analysis
  N = 350,                  # sample size
  r.x1x2.y = .15,           # interaction effect to test (correlation between x1*x2 and y)
  r.x1.y = .2,              # correlation between x1 and y
  r.x2.y = .1,              # correlation between x2 and y
  r.x1.x2 = .2              # correlation between x1 and x2
)
#> [1] "Checking for errors in inputs..."

test_power
#>         pwr shape
#> 1 0.8055776  0.75
```

We see that we have 80% power.

We can also use simulations to estimate power. Simulations are
particularly useful because they can account for non-normal data,
including variable skew, binary variables, and likert variables.

**NB** In all these examples we use 1000 simulations for speed
(`n.iter = 1000`), but for robust results we recommend 10,000
simulations (`n.iter = 10000`).

``` r
test_power<-power_interaction(
  n.iter = 1000,            # number of simulations per unique combination of input parameters
  alpha = 0.05,             # alpha, for the power analysis
  N = 350,                  # sample size
  r.x1x2.y = .15,           # interaction effect to test (correlation between x1*x2 and y)
  r.x1.y = .2,              # correlation between x1 and y
  r.x2.y = .1,              # correlation between x2 and y  
  r.x1.x2 = .2,             # correlation between x1 and x2 
  k.x1 =  2,                # x1 is binary (has 2 levels)
  seed = 581827             # seed, for reproducibility 
)
#> [1] "Checking for errors in inputs..."
#> [1] "Adjusting correlations for variable transformations..."
#> [1] "Performing 1000 simulations"

test_power
#>     N   pwr
#> 1 350 0.834
```

The simulation estimates 83% power - it’s accuracy will increase with
more iterations.

## Changelog

-   **0.1.0.5** *2/17/2022* Adds functions to simulate ordinal variables
    (e.g., a likert scale). Also, a speed-up and minor bug-fixes.
-   **0.1.0.6** *6/22/2022* Added function for analytic power.
