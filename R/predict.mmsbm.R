#' Predict edges based on estimated mmsbm model
#'
#' The function produces expected posterior edges based  
#' on estimated parameters and (optionally new) predictor data 
#'
#' @param fm Object of class \code{mmsbm}.
#' @param new.data.dyad An optional \code{data.frame} object. 
#' @param new.data.monad An optional \code{data.frame} object. 
#' @param forecast Boolean. Should prediction forcast one step into the future? Defaults to FALSE.
#' @param type Character string. The default is to use the linear predictor of edges. The alternative
#'     "response" returns predicted probabilities.    
#'     
#' @return If \code{new.data.dyad = NULL}, vector of length \code{nrow(fm$dyadic.data)}. Else, vector of length \code{nrow(new.data.dyad)}.
#'  
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' 
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
    sid <- fm$forms$senderID
    rid <- fm$forms$receiverID
    dyad <- new.data.dyad
    if(is.null(fm$forms$timeID)){
      tid <- "(tid)"
      dyad[,tid] <- 1
    } else {
      tid <- fm$forms$timeID
    }
  } else {
    sid <- "(sid)"
    rid <- "(rid)"
    tid <- "(tid)"
    dyad <- fm$dyadic.data
  }
  dform <- fm$forms$formula.dyad
  if(any(grepl("missing", names(fm$DyadCoef)))){
    dform <- update(as.formula(fm$forms$formula.dyad), 
                    paste(c("~ .", names(fm$DyadCoef)[grep("missing", 
                                    names(fm$DyadCoef))]), collapse=" + "))
  }
  X_d <- model.matrix(eval(dform), dyad)[, -1]
  if(length(fm$DyadCoef)==0){
    fm$DyadCoef <- as.vector(0)
  }
  if(!is.null(new.data.monad)){
    nid <- ifelse(fm$forms$nodeID %in% colnames(new.data.monad), fm$forms$nodeID, "(nid)")
    monad <- new.data.monad
    if(is.null(fm$forms$timeID)){
      tid <- "(tid)"
      monad[,tid] <- 1
    } else {
      tid <- fm$forms$timeID
    }
  } else {
    nid <- "(nid)"
    tid <- "(tid)"
    monad <- fm$monadic.data
  }
  n_blk <- fm$n_blocks
  mform <- fm$forms$formula.monad
  if(any(grepl("missing", rownames(fm$MonadCoef)))){
    mform <- update(as.formula(fm$forms$formula.monad), 
                    paste(c("~ .", rownames(fm$MonadCoef)[grep("missing", 
                                      rownames(fm$MonadCoef))]), collapse=" + "))
  }
  if(!parametric_mm){
    p <- fm$MixedMembership
  } else {
    if(is.null(fm$forms$formula.monad)){
      X_m <- model.matrix(~ 1, data = monad)
    } else {
      X_m <- model.matrix(eval(mform), monad)
    }
    alpha <- .pi.hat(X_m, fm$MonadCoef)
    if(forecast){
      ts <- unique(monad[,tid])
      new_kappa <- as.matrix(fm$Kappa[,ncol(fm$Kappa)] %*% .mpower(fm$TransitionKernel, forecast))
      new_kappa1 <- matrix(new_kappa, nrow=ncol(new_kappa), ncol=nrow(monad[monad[,tid]==ts[1],]),byrow=FALSE)
      if(length(ts) > 1){
        for(t in 2:length(ts)){
          new_kappa <- rbind(new_kappa, new_kappa[t-1,] %*% .mpower(fm$TransitionKernel, forecast))
          new_kappa1 <- cbind(new_kappa1, matrix(new_kappa[t,], nrow=ncol(new_kappa), ncol=nrow(monad[monad[,tid]==ts[t],]),byrow=FALSE))
        }
      }
      p <- .e.pi(lapply(alpha, function(x)prop.table(x, 2)),
                 new_kappa1)
    } else {
      if(!(tid %in% colnames(monad))){tid <- "(tid)"}
      p <- .e.pi(lapply(alpha, function(x)prop.table(x, 2)),
                 fm$Kappa[,as.character(monad[,tid])])
    }
  }
  
  s_ind <- match(paste(dyad[,sid],dyad[,tid],sep="@"), 
                 paste(monad[,nid],monad[,tid],sep="@"))
  r_ind <- match(paste(dyad[,rid],dyad[,tid],sep="@"), 
                 paste(monad[,nid],monad[,tid],sep="@"))
  pi_s <- p[,s_ind]
  pi_r <- p[,r_ind]
  n_dyad <- nrow(X_d)
  eta_dyad <- array(0.0, n_dyad)
  for(i in 1:n_dyad){ 
      eta_dyad[i] <- pi_s[,i] %*% fm$BlockModel %*% pi_r[,i]
  }
  eta_dyad <- eta_dyad + c(X_d %*% (fm$DyadCoef))
  if(type=="link"){
    return(c(eta_dyad))
  } else {
    return(plogis(c(eta_dyad)))
  }
}

