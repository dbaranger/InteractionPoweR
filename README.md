
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InteractionPoweR

<!-- badges: start -->

<!-- badges: end -->

InteractionPoweR is an R package for running power analyses for
interactions between continuous variables (also known as a moderation
analysis). The main function is `power_interaction()`, which performs
the power analysis. This is done via Monte Carlo simulation.
`power_estimate()` helps to interpret the results of the power
simulation, and `plot_power_curve()` and `plot_simple_sploe()` generate
plots to visualize the results. The function `generate_interaction()`
simulates a single data set drawn from the specified population-level
effects and `plot_interaction()` can be used to visualize the simulated
data.

## Installation

You can install InteractionPoweR from github with:

``` r
install.packages("devtools")
devtools::install_github("dbaranger/InteractionPoweR")
```

## Examples

### Example 1

The simplest use-case is when all the input parameters are known. We
know the population-level correlation between our predictors (x1 and x2)
and our outcome, we have a smallest effect size of interest in mind for
our interaction effect size, and our sample size is already set (maybe
we are conducting secondary data analysis). Power can be determined with
a single command:

``` r
library(InteractionPoweR)

test_power<-power_interaction(
  n.iter = 1000,            # number of simulations per unique combination of input parameters
  cl = 6,                   # number of cores for parallel processing (strongly recommended)
  alpha = 0.05,             # alpha, for the power analysis
  N = 350,                  # sample size
  r.x1x2.y = .15,           # interaction effect to test (correlation between x1*x2 and y)
  r.x1.y = .2,              # correlation between x1 and y
  r.x2.y = .1,              # correlation between x2 and y
  r.x1.x2 = .2              # correlation between x1 and x2
)
#> [1] "Performing 1000 tests"
test_power
#>     pwr     min.lwr   min.upr   max.lwr   max.upr
#> 1 0.822 -0.07767766 0.2013109 0.1825631 0.4779925
```

We see that we have \~80% power to detect the effect of interest.
`min.lwr`, `min.upr` etc reflect the range of simple-slopes you would
expect to see.

It can be hard to know what interaction correlations mean in terms of
how the data will look. To help users interpret interaction effects, we
provide a simple interface for simulating single data sets and plotting
them.

``` r
sample_data<-generate_interaction(N=350,r.x1x2.y =.15,r.x1.y = .2, r.x2.y = .1, r.x1.x2 = .2)
plot_interaction(data = sample_data)
```

<img src="man/figures/README-example2-1.png" width="100%" />

### Example 2

In this example, we know the population-level correlation between each
of our predictors (x1 and x2) and our outcome (y), as well as the
correlation between the two predictors. We are interested in
interactions within a certain range, and wish to know what sample size
we would need to detect those interactions.

``` r
library(tictoc)
tic()

test_power<-power_interaction(
  n.iter = 1000,            # number of simulations per unique combination of input parameters
  cl = 6,                   # number of cores for parallel processing (strongly recommended)
  alpha = 0.05,             # alpha, for the power analysis
  N = seq(20,400,by=20),    # range of sample sizes to test
  r.x1x2.y = c(.18,.2,.22), # range of interaction effects to test
  r.x1.y = .2,              # correlation between x1 and y
  r.x2.y = .1,              # correlation between x2 and y
  r.x1.x2 = .2              # correlation between x1 and x2
)
#> [1] "Performing 60000 tests"
toc()
#> 194.25 sec elapsed
```

The results of this analysis can be hard to interpret just by looking at
the output. Instead, we recommend visualizing them using
`plot_power_curve()`. In this case, the horizontal line on the plot will
be at our target power of 90%:

``` r
plot_power_curve(test_power,power_target = .9)
```

<img src="man/figures/README-example4-1.png" width="100%" />

The function `power_estimate()` can be used to estimate where the
power\_curve for each interaction effect size crosses our 90% line:

``` r
power_estimate(test_power,power_target = .9,x = "N")
#>   r.x1x2.y estimate
#> 1     0.18 317.6146
#> 2     0.20 245.1963
#> 3     0.22 208.2832
```

We can see that depending on the specific effect size we hope to detect,
we would need between N\~200 and N\~300 participants.
