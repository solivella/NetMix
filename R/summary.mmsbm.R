#' Summarize 'mmsbm' object
#'
#' The function summarizes the output of a dynMMSBM model object
#'
#' @param object An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param ... Currently ignored
#' @return List with named components:
#'     \describe{
#'       \item{N}{Total number of dyad-time period observations.}
#'       \item{Number of Clusters}{Number of latent groups included in the dynMMSBM model.}
#'       \item{Percent of Observations in Each Cluster}{Average membership in each latent group, across all node-time periods.}
#'       \item{Edge Formation Probabilities}{\code{n.groups} by \code{n.groups} matrix of estimated edge formation probabilities between latent groups.}
#'       \item{Dyadic Coefficients}{Vector of estimated coefficient values for dyadic covariates.}
#'       \item{Monadic Coefficients}{Array of estimated coefficient values for monadic covariates. Has \code{n.groups} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{Markov State Probabilities}{Average HMM state probabilities across all time periods.}
#'     }
#'     
#' @method summary mmsbm
#' 
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 3 groups
#' set.seed(123)
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 3)
#' 
#' ## Summarize estimated model
#' summary(lazega_mmsbm)






summary.mmsbm <- function(object, ...){
  summ <- list(nrow(object$dyadic.data), ncol(object$BlockModel), 
               rowMeans(object$MixedMembership),
               exp(object$BlockModel) / (1 + exp(object$BlockModel)), 
               object$DyadCoef, object$MonadCoef, rowMeans(object$Kappa))
  names(summ) <- c("N", "Number of Clusters", "Percent of Observations in Each Cluster",
                   "Edge Formation Probabilities", "Dyadic Coefficients", "Monadic Coefficients",
                   "Markov State Probabilities")
  print(summ)
}
