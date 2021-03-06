% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.mmsbm.R
\name{simulate.mmsbm}
\alias{simulate.mmsbm}
\title{Simulate a complete sociomatrix from an \code{mmsbm} object}
\usage{
\method{simulate}{mmsbm}(
  object,
  nsim = 1,
  seed = NULL,
  new.data.dyad = NULL,
  new.data.monad = NULL,
  ...
)
}
\arguments{
\item{object}{An object of class \code{mmsbm}, a result of a call to \code{mmsbm}}

\item{nsim}{Number of networks to simulate}

\item{seed}{RNG seed.}

\item{new.data.dyad}{An optional \code{data.frame} object. If not \code{NULL}, use these 
dyadic predictor values instead of those used to fit the original model.}

\item{new.data.monad}{An optional \code{data.frame} object. See \code{new.data.dyad}.}

\item{...}{Currently ignored}
}
\value{
List of length \code{nsim} of simulated networks. 
        If \code{new.data.dyad = NULL}, each element is a vector of length \code{nrow(object$dyadic.data)}. 
        Else, vector of length \code{nrow(new.data.dyad)}. If \code{seed} is not NULL, return object
        includes its value as attribute "seed".
}
\description{
The function generates one sample network from the posterior predictive of the model represented by a fitted \code{mmsbm} object.
}
\examples{
library(NetMix)
## Load datasets
data("lazega_dyadic")
data("lazega_monadic")
## Estimate model with 2 groups
lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
                      ~  School + Practice + Status,
                      senderID = "Lawyer1",
                      receiverID = "Lawyer2",
                      nodeID = "Lawyer",
                      data.dyad = lazega_dyadic,
                      data.monad = lazega_monadic,
                      n.blocks = 2,
                      mmsbm.control = list(seed = 123, 
                                           conv_tol = 1e-2,
                                           hessian = FALSE))

## Simulate 5 new networks
lazega_sim <- simulate(lazega_mmsbm, nsim = 5, seed = 123)

}
\author{
Santiago Olivella (olivella@unc.edu), Adeline Lo (aylo@wisc.edu), Tyler Pratt (tyler.pratt@yale.edu), Kosuke Imai (imai@harvard.edu)
}
