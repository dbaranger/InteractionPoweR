% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformation_functions.R
\name{compute_adjustment}
\alias{compute_adjustment}
\title{compute_adjustment}
\usage{
compute_adjustment(
  r.x1.y,
  r.x2.y,
  r.x1x2.y,
  r.x1.x2,
  N.adjustment = 1e+06,
  tol = 0.005,
  iter = 10,
  skew.x1,
  skew.x2,
  skew.y,
  transform.x1,
  transform.x2,
  transform.y
)
}
\arguments{
\item{r.x1.y}{Internal use only}

\item{r.x2.y}{Internal use only}

\item{r.x1x2.y}{Internal use only}

\item{r.x1.x2}{Internal use only}

\item{N.adjustment}{Internal use only}

\item{tol}{Internal use only}

\item{iter}{Internal use only}

\item{skew.x1}{Internal use only}

\item{skew.x2}{Internal use only}

\item{skew.y}{Internal use only}

\item{transform.x1}{Internal use only}

\item{transform.x2}{Internal use only}

\item{transform.y}{Internal use only}
}
\value{
Correlation adjustments.
}
\description{
Computes how much variable correlations need to be adjusted so that they have the desired correlation structure after transformation. Intended for internal use only.
}
\examples{
\dontrun{
compute_adjustment(r.x1.y = .2,r.x2.y = .2,r.x1x2.y = .1,r.x1.x2 = .2,transform.y = "binary")
}
}