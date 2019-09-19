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
  summ <- list(nrow(object$dyadic.data), ncol(object$BlockModel), 
               rowMeans(object$MixedMembership),
               exp(object$BlockModel) / (1 + exp(object$BlockModel)), 
               object$DyadCoef, object$MonadCoef, rowMeans(object$Kappa))
  names(summ) <- c("N", "Number of Clusters", "Percent of Observations in Each Cluster",
                   "Edge Formation Probabilities", "Dyadic Coefficients", "Monadic Coefficients",
                   "Markov State Probabilities")
  if("vcov_dyad" %in% names(object)){
    summ$`Dyadic Coefficients` <- cbind(object$DyadCoef,
                                        sqrt(diag(object$vcov_dyad)))
    colnames(summ$`Dyadic Coefficients`) <- c("Coefficient", "Std. Error")
    
    colnames(summ$`Monadic Coefficients`) <- paste(colnames(summ$`Monadic Coefficients`), "Coefficient")
    mse <- sqrt(diag(object$vcov_monad))
    z <- array(dim=c(nrow(object$MonadCoef), ncol(object$MonadCoef)*2, object$n_states))
    dimnames(z)[c(1,3)] <- dimnames(summ$`Monadic Coefficients`)[c(1,3)]
    for(i in 1:object$n_states){
      m <- cbind(summ$`Monadic Coefficients`[,,i],
                 matrix(nrow=nrow(object$MonadCoef), ncol=object$n_blocks,
                        mse[grep(paste("State", i), names(mse))], 
                        dimnames = list(NULL, paste(colnames(object$MonadCoef), "Std. Error"))))
      z[,,i] <- m[,c(rbind(1:ncol(object$MonadCoef), (1:ncol(object$MonadCoef))+ncol(object$MonadCoef)))]
      colnames(z) <- colnames(m[,c(rbind(1:ncol(object$MonadCoef), (1:ncol(object$MonadCoef))+ncol(object$MonadCoef)))])
    }
    summ$`Monadic Coefficients` <- z
  }
  print(summ)
}
