% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_functions.R
\name{get_psi}
\alias{get_psi}
\title{get_psi}
\usage{
get_psi(df, X, Y, Astar, estimand, sl.lib)
}
\arguments{
\item{df}{A dataframe containing the outcome, treatment and covariates}

\item{X}{A vector of covariates}

\item{Y}{A vector of outcomes}

\item{estimand}{The estimand of interest (e.g. ATE, ATT, ATC)}

\item{sl.lib}{A vector of super learner libraries to use}

\item{A}{A vector of treatments}
}
\value{
A list containing the AIPW estimate of the control variates
}
\description{
get_psi
}
