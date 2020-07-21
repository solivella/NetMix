#' Summarize fitted mmsbmB
#'
#' The function summarizes the output of a mmsbmB model object
#'
#' @param fm An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}.
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

summary.mmsbmB <- function(object,...){
  summ <- list("Number of Dyads" = nrow(object$dyadic.data),
               "Number of Family 1 Blocks" = nrow(object$BlockModel), 
               "Number of Family 2 Blocks" = ncol(object$BlockModel), 
               "Percent of Observations in Each Family 1 Block" = rowMeans(object$`MixedMembership 1`),
               "Percent of Observations in Each Family 2 Block" = rowMeans(object$`MixedMembership 2`),
               "Blockmodel Matrix" = exp(t(object$BlockModel)) / (1 + exp(t(object$BlockModel))), ##### AFTER FIX OF return t(B) change this back
               "Family 1 Monadic Coefficients" = object$MonadCoef1,
               "Family 2 Monadic Coefficients" = object$MonadCoef2)  
  
  if(nrow(summ$`Family 1 Monadic Coefficients`)>1){
    rownames(summ$`Family 1 Monadic Coefficients`)<-colnames(model.matrix(object$forms$formula.monad1,object$monadic1.data))
  }else{rownames(summ$`Family 1 Monadic Coefficients`)<-"Intercept"}
  
  if(nrow(summ$`Family 2 Monadic Coefficients`)>1){
    rownames(summ$`Family 2 Monadic Coefficients`)<-colnames(model.matrix(object$forms$formula.monad2,object$monadic2.data))
  }else{rownames(summ$`Family 2 Monadic Coefficients`)<-"Intercept"}
  
  if(length(object$DyadCoef)){
    summ$`Dyadic Coefficients` <- object$DyadCoef
    rownames(summ$`Dyadic Coefficients`) <- colnames(model.matrix(object$forms$formula.dyad,object$dyadic.data))[-1]
  }
  if(object$n_states > 1){
    summ$`Markov State Probabilities` <-  rowMeans(object$Kappa)
  }
  
  # Currently edit out below, no hessian in bip set up yet
  # if(object$forms$hessian){
  #   if("vcov_dyad" %in% names(object)){
  #     summ$`Dyadic Coefficients` <- cbind(object$DyadCoef,
  #                                         sqrt(diag(object$vcov_dyad)))
  #     colnames(summ$`Dyadic Coefficients`) <- c("Coefficient", "Std. Error")
  #   }
  #   
  #  mse <- sqrt(diag(object$vcov_monad))
  # summ$`Monadic Coefficients` <- cbind(c(summ$`Monadic Coefficients`),
  #                                      mse)
  # colnames(summ$`Monadic Coefficients`) <- c("Coefficient", "Std. Error")
  # rownames(summ$`Monadic Coefficients`) <- rownames(object$vcov_monad)
  
  #}
  print(summ)
}
