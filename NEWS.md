# InteractionPoweR 0.2.2

-   Adds power analysis function for 3-way interactions -
    `power_interaction_3way_r2()`.
-   Adds convenience function for visualizing the correlation matrix
    from a 3-way interaction - `cor.mat.3way()`.
-   Adds convenience function for visualizing the simple slopes from a
    3-way interaction - `simple.slopes.3way()`.
-   Adds power analysis function for 2-way interactions with
    covariates - `generate.interaction.cov.input()` and
    `power_interaction_r2_covs()`.
-   Minor bug fixes.
-   Removes dependency on the {MASS} package.

# InteractionPoweR 0.2.1

-   Fix bug in `power_interaction_r2()` which effected analyses with a
    high ‘r.x1.x2’ correlation.
-   Fix warning in `plot_interaction()`.
-   Adds additional descriptive error messages.

# InteractionPoweR 0.2.0

-   Removes options for specifying skew, which contained a bug.
-   Also a speed-up for the simulations.

# InteractionPoweR 0.1.1

-   First CRAN release.

# InteractionPoweR 0.1.0.6

-   Added `power_interaction_r2()` function for analytic power.

# InteractionPoweR 0.1.0.5

-   Adds options to simulate ordinal variables (e.g., a likert scale) in
    `power_interaction()` using the ‘k.x1’, ‘k.x2’, and ‘k.y’ options.
-   Also a speed-up and minor bug-fixes.
