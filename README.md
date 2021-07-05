
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InteractionPoweR

<!-- badges: start -->
<!-- badges: end -->

InteractionPoweR is an R package for running power analyses for
interactions in cross-sectional data sets between continuous and/or
binary variables (also known as a moderation analysis). The main
function is `power_interaction()`, which performs the power analysis.
This is done via Monte Carlo simulation. `power_estimate()` helps to
interpret the results of the power simulation, and `plot_power_curve()`
and `plot_simple_slope()` generate plots to visualize the results. The
function `generate_interaction()` simulates a single data set drawn from
the specified population-level effects and `plot_interaction()` can be
used to visualize the simulated data.

**See the [website](https://dbaranger.github.io/InteractionPoweR/) for
further documentation and examples.**

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
#> 18.92 sec elapsed
test_power
#> # A tibble: 1 x 1
#>     pwr
#>   <dbl>
#> 1 0.809
```

We see that we have 80.9% power to detect the effect of interest.
