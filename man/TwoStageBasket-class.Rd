% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\docType{class}
\name{TwoStageBasket-class}
\alias{TwoStageBasket-class}
\alias{TwoStageBasket}
\title{Class TwoStageBasket}
\description{
OneStageBasket is an S4 class. An object of this class contains the most
important design features of a single-stage basket trial.
}
\details{
This class implements a two-stage basket trial based on the design
proposed by Fujikawa et al.
}
\section{Slots}{

\describe{
\item{\code{k}}{The number of baskets.}

\item{\code{shape1}}{First common shape parameter of the beta prior.}

\item{\code{shape2}}{Second common shape parameter of the beta prior.}

\item{\code{theta0}}{A common probability under the null hypothesis.}
}}

