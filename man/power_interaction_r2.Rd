% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_interaction_r2.R
\name{power_interaction_r2}
\alias{power_interaction_r2}
\title{Analytic power analysis for interactions}
\usage{
power_interaction_r2(
  N,
  r.x1.y,
  r.x2.y,
  r.x1x2.y,
  r.x1.x2,
  rel.x1 = 1,
  rel.x2 = 1,
  rel.y = 1,
  alpha = 0.05,
  detailed_results = FALSE
)
}
\arguments{
\item{N}{Sample size. Must be a positive integer. Has no default value. Can be a single value or a vector of values.}

\item{r.x1.y}{Pearson's correlation between x1 and y. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.}

\item{r.x2.y}{Pearson's correlation between x2 and y. Must be between -1 and 1. Assumed to be the 'moderator' in some functions. Has no default value. Can be a single value or a vector of values.}

\item{r.x1x2.y}{Pearson's correlation between the interaction term x1x2 (x1 * x2) and y. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.}

\item{r.x1.x2}{Pearson's correlation between x1 and x2. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.}

\item{rel.x1}{Reliability of x1 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.}

\item{rel.x2}{Reliability of x2 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.}

\item{rel.y}{Reliability of xy (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.}

\item{alpha}{The alpha. At what p-value is the interaction deemed significant? Default is 0.05.}

\item{detailed_results}{Default is FALSE. Should detailed results be reported?}
}
\value{
A data frame containing the power for each unique setting combination.
}
\description{
Power analysis for interaction models, computed via change in R2. Valid for interactions with continuous, normally distributed, variables.
}
\examples{
power_interaction_r2(N=seq(100,300,by=10),r.x1.y=0.2, r.x2.y=.2,r.x1x2.y=0.2,r.x1.x2=.2)

}
