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
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (adelinel@@princeton.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 2 groups
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(seed = 123,
#'                                            hessian = TRUE))
#' 
#' ## Summarize estimated model
#' summary(lazega_mmsbm)
#' 

summary.mmsbm <- function(object, ...){
  summ <- list("Number of Dyads" = nrow(object$dyadic.data),
               "Number of Blocks" = ncol(object$BlockModel), 
               "Percent of Observations in Each Block" = rowMeans(object$MixedMembership),
               "Blockmodel Matrix" = exp(object$BlockModel) / (1 + exp(object$BlockModel)),
               "Monadic Coefficients" = object$MonadCoef)
  if(length(object$DyadCoef)){
    summ$`Dyadic Coefficients` <- object$DyadCoef
  }
  if(object$n_states > 1){
    summ$`Markov State Probabilities` <-  rowMeans(object$Kappa)
  }

  if(object$forms$hessian){
    if("vcov_dyad" %in% names(object)){
      summ$`Dyadic Coefficients` <- cbind(object$DyadCoef,
                                          sqrt(diag(object$vcov_dyad)))
      colnames(summ$`Dyadic Coefficients`) <- c("Coefficient", "Std. Error")
    }
    
    mse <- sqrt(diag(object$vcov_monad))
    summ$`Monadic Coefficients` <- cbind(c(summ$`Monadic Coefficients`),
                                         mse)
    colnames(summ$`Monadic Coefficients`) <- c("Coefficient", "Std. Error")
    rownames(summ$`Monadic Coefficients`) <- rownames(object$vcov_monad)

  }
  print(summ)
}
