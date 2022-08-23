## Resubmission

This is a resubmission. In this version I have addressed the recommended
changes to the package:

-   Added () behind function names in the DESCRIPTION.
-   Reformatted the reference in the DESCRIPTION.
-   Use TRUE or FALSE instead of T/F.
-   Removed or replaced with .
-   Removed information messages that cannot be suppressed (replaced
    with message() or warning()).
-   No longer sets a seed within functions.

## R CMD check results

There were no ERRORs or WARNINGs. NOTEs are related to the spelling of
author’s names.

## Downstream dependencies

There are no downstream dependencies.

## Test environments

Passes R CMD check with no ERRORs, WARNINGs, or NOTEs on:

1.  Windows Server 2022, R-release, 32/64 bit
2.  Windows Server 2022, R-devel, 64 bit
3.  Apple Silicon (M1), macOS 11.6 Big Sur, R-release
4.  macOS 10.13.6 High Sierra, R-release, brew
