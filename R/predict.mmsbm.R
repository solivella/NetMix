#' Predict edges based on estimated mmsbm model
#'
#' The function produces expected posterior edges based  
#' on estimated parameters and (optionally new) predictor data 
#'
#' @param object Object of class \code{mmsbm}.
#' @param new.data.dyad An optional \code{data.frame} object. 
#' @param new.data.monad An optional \code{data.frame} object. 
#' @param forecast Boolean. Should prediction forcast one step into the future? Defaults to FALSE.
#' @param type Character string. The default is to use the linear predictor of edges. The alternative
#'     "response" returns predicted probabilities. The alternative "mm" returns predicted mixed-membership vectors.  
#' @param ... Currently ignored  
#'     
#' @return If \code{new.data.dyad = NULL}, vector of length \code{nrow(object$dyadic.data)}. Else, vector of length \code{nrow(new.data.dyad)}.
#' 
#' @method predict mmsbm
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
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(seed = 123,
#'                                            conv_tol = 1e-2, 
#'                                            hessian = FALSE))
#' 
#' ## Get in-sample predicted edge probabilities
#' lazega_preds <- predict(lazega_mmsbm, type = "response")
#' 



predict.mmsbm <- function(object, 
                          new.data.dyad = NULL,
                          new.data.monad  = NULL, 
                          forecast = FALSE,
                          type = c("link", "response","mm"),
                          ...)
{
  type <- match.arg(type)
  if(!is.null(new.data.dyad)){
    sid <- object$forms$senderID
    rid <- object$forms$receiverID
    dyad <- new.data.dyad
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      dyad[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    sid <- "(sid)"
    rid <- "(rid)"
    tid <- "(tid)"
    dyad <- object$dyadic.data
  }
  dform <- object$forms$formula.dyad
  if(any(grepl("missing", names(object$DyadCoef)))){
    dform <- update(as.formula(dform), 
                    paste(c("~ .", names(object$DyadCoef)[grep("missing", 
                                    names(object$DyadCoef))]), collapse=" + "))
  }
  X_d <- model.matrix(eval(dform), dyad)
  if(length(object$DyadCoef)==0){
    object$DyadCoef <- as.vector(0)
  } else {
    object$DyadCoef <- c(0, object$DyadCoef)
  }
  if(!is.null(new.data.monad)){
    nid <- ifelse(object$forms$nodeID %in% colnames(new.data.monad), object$forms$nodeID, "(nid)")
    monad <- new.data.monad
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      monad[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
    C_mat <- matrix(0, ncol = object$forms$n.blocks, nrow = nrow(monad))
  } else {
    nid <- "(nid)"
    tid <- "(tid)"
    monad <- object$monadic.data
    C_mat <- object$CountMatrix
  }
  n_blk <- object$n_blocks
  mform <- object$forms$formula.monad
  if(any(grepl("missing", rownames(object$MonadCoef)))){
    mform <- update(as.formula(mform), 
                    paste(c("~ .", rownames(object$MonadCoef)[grep("missing", 
                                      rownames(object$MonadCoef))]), collapse=" + "))
  }

    if(is.null(mform)){
      X_m <- model.matrix(~ 1, data = monad)
    } else {
      X_m <- model.matrix(eval(mform), monad)
    }
    alpha <- .compute.alpha(X_m, object$MonadCoef)
    if(forecast){
      ts <- unique(monad[,tid])
      new_kappa <- as.matrix(object$Kappa[,ncol(object$Kappa)] %*% .mpower(object$TransitionKernel, forecast))
      new_kappa1 <- matrix(new_kappa, nrow=ncol(new_kappa), ncol=nrow(monad[monad[,tid]==ts[1],]),byrow=FALSE)
      if(length(ts) > 1){
        for(t in 2:length(ts)){
          new_kappa <- rbind(new_kappa, new_kappa[t-1,] %*% .mpower(object$TransitionKernel, forecast))
          new_kappa1 <- cbind(new_kappa1, matrix(new_kappa[t,], nrow=ncol(new_kappa), ncol=nrow(monad[monad[,tid]==ts[t],]),byrow=FALSE))
        }
      }
      p <- .e.pi(alpha, new_kappa1, C_mat)
    } else {
      if(!(tid %in% colnames(monad))){tid <- "(tid)"}
      p <- .e.pi(alpha, object$Kappa[,as.character(monad[,tid])], C_mat)
    }
  if(type == "mm"){
    return(p)
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
      eta_dyad[i] <- pi_s[,i] %*% object$BlockModel %*% pi_r[,i]
  }
  eta_dyad <- eta_dyad + c(X_d %*% (object$DyadCoef))
  if(type == "link"){
    return(c(eta_dyad))
  } else if (type == "response") {
    return(plogis(c(eta_dyad)))
  }
}

