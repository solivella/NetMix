#' Summarize fitted mmsbmB
#'
#' The function summarizes the output of a mmsbmB model object
#'
#' @param object An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}.
#' @return List with named components:
#'     \describe{
#'       \item{N}{Total number of dyad-time period observations.}
#'       \item{Number of Clusters 1}{Number of latent groups included in the mmsbmB model for Family 1.}
#'       \item{Number of Clusters 2}{Number of latent groups included in the mmsbmB model for Family 2.}
#'       \item{Percent of Observations in Each Cluster 1}{Average membership in each latent group, across all node-time periods, Family 1.}
#'       \item{Percent of Observations in Each Cluster 2}{Average membership in each latent group, across all node-time periods, Family 2.}
#'       \item{Edge Formation Probabilities}{\code{n.groups1} by \code{n.groups2} matrix of estimated edge formation probabilities between latent groups.}
#'       \item{Dyadic Coefficients}{Vector of estimated coefficient values for dyadic covariates.}
#'       \item{Monadic Coefficients 1}{Array of estimated coefficient values for monadic covariates in Family 1. Has \code{n.groups1} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{Monadic Coefficients 2}{Array of estimated coefficient values for monadic covariates in Family 2. Has \code{n.groups2} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{Markov State Probabilities}{Average HMM state probabilities across all time periods. Not currently in usage.}
#'     }
#'     
#' @method summary mmsbmB
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)


summary.mmsbmB <- function(object,...){
  summ <- list("Number of Dyads" = nrow(object$dyadic.data),
               "Number of Family 1 Blocks" = nrow(object$BlockModel), 
               "Number of Family 2 Blocks" = ncol(object$BlockModel), 
               "Percent of Observations in Each Family 1 Block" = rowMeans(object$MixedMembership1),
               "Percent of Observations in Each Family 2 Block" = rowMeans(object$MixedMembership2),
               "Blockmodel Matrix" = exp(object$BlockModel) / (1 + exp(object$BlockModel)), 
               "Monadic Coefficients 1" = object$MonadCoef1,
               "Monadic Coefficients 2" = object$MonadCoef2)  
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
    #family 1
    mse1 <- sqrt(diag(object$vcov_monad1))
    summ$`Monadic Coefficients 1` <- cbind(c(summ$`Monadic Coefficients 1`),
                                         mse1)
    colnames(summ$`Monadic Coefficients 1`) <- c("Coefficient", "Std. Error")
    rownames(summ$`Monadic Coefficients 1`) <- rownames(object$vcov_monad1)
    #family 2
    mse2 <- sqrt(diag(object$vcov_monad2))
    summ$`Monadic Coefficients 2` <- cbind(c(summ$`Monadic Coefficients 2`),mse2)
    colnames(summ$`Monadic Coefficients 2`) <- c("Coefficient", "Std. Error")
    rownames(summ$`Monadic Coefficients 2`) <- rownames(object$vcov_monad2)
    
  }
  
  
  print(summ)
}
