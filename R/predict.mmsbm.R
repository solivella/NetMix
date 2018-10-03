#' Predict edges based on estimated mmsbm model
#'
#' The function produces expected posterior edges based  
#' on estimated parameters and (optionally new) predictor data 
#'
#' @param fm Object of class \code{mmsbm}.
#' @param new.data.dyad An optional \code{data.frame} object. 
#' @param new.data.monad An optional \code{data.frame} object. 
#' @param forecast Boolean. Defaults to FALSE.
#' @param type Character string. The default is to use the linear predictor of edges. The alternative
#'     "response" returns predicted probabilities.    
#'     
#' @return If \code{new.data.dyad = NULL}, vector of length \code{nrow(fm$dyadic.data)}. Else, vector of length \code{nrow(new.data.dyad)}.
#'  
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' @example tests/Examples/MIDColdWar.R
#' 
predict.mmsbm <- function(fm, 
                          new.data.dyad = NULL,
                          new.data.monad  = NULL, 
                          parametric_mm = FALSE,
                          forecast = FALSE,
                          type = c("link", "response"))
{
  type <- match.arg(type)
  if(forecast & !parametric_mm){
    stop("Must use parametric mixed-membership terms when forecasting.")
  }
  if(!is.null(new.data.dyad)){
    sid <- fm$call$senderID
    rid <- fm$call$receiverID
    dyad <- new.data.dyad
    if(is.null(fm$call$timeID)){
      tid <- "(tid)"
      dyad[,tid] <- 1
    } else {
      tid <- fm$call$timeID
    }
  } else {
    sid <- "(sid)"
    rid <- "(rid)"
    tid <- "(tid)"
    dyad <- fm$dyadic.data
  }
  X_d <- model.matrix(eval(fm$call$formula.dyad), dyad)[, -1]
  if(length(fm$DyadCoef)==0){
    fm$DyadCoef <- as.vector(0)
  }
  if(!is.null(new.data.monad)){
    nid <- ifelse(fm$call$nodeID %in% colnames(new.data.monad), fm$call$nodeID, "(nid)")
    monad <- new.data.monad
    if(is.null(fm$call$timeID)){
      tid <- "(tid)"
      monad[,tid] <- 1
    } else {
      tid <- fm$call$timeID
    }
  } else {
    nid <- "(nid)"
    tid <- "(tid)"
    monad <- fm$monadic.data
  }
  n_blk <- fm$n_blocks
  if(!parametric_mm){
    p <- fm$MixedMembership
  } else {
    if(is.null(fm$call$formula.monad)){
      X_m <- model.matrix(~ 1, data = monad)
    } else {
      X_m <- model.matrix(eval(fm$call$formula.monad), monad)
    }
    ph <- .pi.hat(X_m, fm$MonadCoef)
    if(forecast){
      new_kappa <- fm$Kappa[,ncol(fm$Kappa)] %*% .mpower(fm$TransitionKernel, forecast)
      p <- .e.pi(ph, new_kappa)
    } else {
      p <- .e.pi(ph, fm$Kappa[,monad[,tid]])
    }
  }
  eta_dyad <- X_d %*% as.matrix(fm$DyadCoef)
  s_ind <- match(dyad[,sid], monad[,nid])
  r_ind <- match(dyad[,rid], monad[,nid])
  for(a in 1:n_blk){ 
    for(b in 1:n_blk){
      eta_dyad <- eta_dyad + (p[a,s_ind]*p[b,r_ind]*fm$BlockModel[a,b])
    }
  }
  if(type=="link"){
    return(c(eta_dyad))
  } else {
    return(plogis(c(eta_dyad)))
  }
}

