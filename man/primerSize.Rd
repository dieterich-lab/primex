% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/settings.R
\name{primerSize}
\alias{primerSize}
\title{Set primer size}
\usage{
primerSize(settings, min = NULL, optimal = NULL, max = NULL)
}
\arguments{
\item{settings}{a settings list created by `p3Settings`}

\item{min}{minimal length (PRIMER_MIN_SIZE)}

\item{optimal}{optimal length (PRIMER_OPT_SIZE)}

\item{max}{maximal length (PRIMER_MAX_SIZE)}
}
\value{
a settings list with updated values
}
\description{
Set primer size
}
\examples{
p3Settings() \%>\%
  primerSize(min = 18, optimal = 20, max = 22) \%>\%
  str()
}
