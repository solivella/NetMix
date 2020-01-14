#' Generate estimated monadic covariate effects for estimated mmsbm model
#'
#' The function estimates the effect of a shift in monadic covariate values on the probability of edge formation in the network. 
#'
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}. 
#' @param cov Character string identifying the monadic covariate to be shifted.  
#' @param shift Numeric value specifying the desired increase or decrease in the monadic covariate.  The monadic predictor will be shifted by this value for all nodes and time periods.
#' @param max.val An optional numeric value specifying the maximum possible value for the monadic covariate.
#'
#'     
#' @return List with named components:
#'     \describe{
#'       \item{Overall Avg. Effect}{Overall average effect of the covariate shift on the predicted probability of edge formation.}
#'       \item{Avg. Effect by Time}{Vector of average effects of the covariate shift on the predicted probability of edge formation for each time period.}
#'       \item{Avg. Effect by Node}{Vector of average effects of the covariate shift on the predicted probability of edge formation for each node.}
#'       \item{Avg. Effect by Dyad}{Vector of average effects of the covariate shift on the predicted probability of edge formation for each node dyad.}
#'       \item{Avg. Effect Dyad-Time}{Vector of estimated effects of the covariate shift on the predicted probability of edge formation for each node dyad-time unit.}
#'     }
#' 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 2 groups
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  Age,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(seed = 123, 
#'                                            hessian = FALSE))
#' 
#' ## Compute effect of decreasing every lawyers' age by 10 years
#' fx_list <- covFX(lazega_mmsbm, cov = "Age", shift = -10)
#' fx_list[["Overall Avg. Effect of Age"]]
#' 



covFX <- function(fm, cov, shift, max.val=FALSE){
  predict.ties <- predict(fm, parametric_mm = TRUE, type="response")
  monadic.data2 <- fm$monadic.data
  monadic.data2[,cov] <- fm$monadic.data[,cov] + shift
  if(!isFALSE(max.val)){
    monadic.data2[which(fm$monadic.data[,cov] == max(fm$monadic.data[,cov])),cov] <- max.val
  }
  predict.ties2 <- predict(fm, new.data.monad=monadic.data2, parametric_mm = TRUE, type="response")
  FX <- list(mean(predict.ties2 - predict.ties), #avg
             tapply(predict.ties2-predict.ties, fm$dyadic.data[,"(tid)"], mean), #time
             sapply(unique(fm$monadic.data[,"(nid)"]), function(x){ #node
               mean((predict.ties2-predict.ties)[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x])}),
             tapply(predict.ties2-predict.ties, paste(fm$dyadic.data[,"(sid)"], fm$dyadic.data[,"(rid)"], sep="_"), mean),#dyad
             predict.ties2 - predict.ties) #dyad-time
  names(FX[[3]]) <- unique(fm$monadic.data[,"(nid)"])
  names(FX[[5]]) <- paste(fm$dyadic.data[,"(sid)"], fm$dyadic.data[,"(rid)"], sep="_")
  names(FX) <- c(paste("Overall Avg. Effect of", cov), paste("Avg. Effect of", cov, "by Time"),
                 paste("Avg. Effect of", cov, "by Node"), paste("Avg. Effect of", cov, "by Dyad"),
                 paste("Effect of", cov, "by Dyad-Time"))
  return(FX)
}
