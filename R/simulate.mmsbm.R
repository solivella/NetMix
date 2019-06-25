#' Simulate a complete sociomatrix from an \code{mmsbm} object
#'
#' The function generates one sample network from the posterior predictive of the model represented by a fitted \code{mmsbm} object.
#'
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param new.data.dyad An optional \code{data.frame} object. If not \code{NULL}, use these 
#'                      dyadic predictor values instead of those used to fit the original model.
#' @param new.data.monad An optional \code{data.frame} object. See \code{new.data.dyad}. 
#' @param parametric_mm Boolean. Should the variational posterior be used for sampling the mixed-memberships (\code{FALSE}), 
#'                      or should the mixed-meberships be formed using the parameters in the monadic regression equation (\code{TRUE})?
#'                      Defaults to \code{FALSE}. If \code{is.null(new.data.monad)=FALSE}, setting this to \code{FALSE} will produce an error.  
#'
#' @return If \code{new.data.dyad = NULL}, vector of length \code{nrow(fm$dyadic.data)}. Else, vector of length \code{nrow(new.data.dyad)}.
#'
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)

simulate.mmsbm <- function(fm, 
                           new.data.dyad = NULL,
                           new.data.monad  = NULL, 
                           parametric_mm = FALSE)
{
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
  if(!is.null(new.data.monad)){
    nid <- fm$forms$nodeID
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
  if(! (tid %in% c(names(dyad), names(monad)))){
    stop("Dynamic model estimated, but no timeID provided in new data.")
  }
  X_d <- model.matrix(fm$forms$formula.dyad, dyad)
  if(is.null(fm$forms$formula.monad)){
    X_m <- model.matrix(~ 1, data = monad)
  } else {
    X_m <- model.matrix(fm$forms$formula.monad, monad)
  }
  if(length(fm$DyadCoef)==0){
    fm$DyadCoef <- as.vector(0)
  } else {
    fm$DyadCoef <- c(0, fm$DyadCoef)
  }
  s_ind <- match(paste(dyad[,sid],dyad[,tid],sep="@"), 
                 paste(monad[,nid],monad[,tid],sep="@"))
  r_ind <- match(paste(dyad[,rid],dyad[,tid],sep="@"), 
                 paste(monad[,nid],monad[,tid],sep="@"))
  n_blk <- fm$n_blocks
  n_dyad <- nrow(X_d)
  
  unique_t <- unique(monad[,tid])
  states <- as.matrix(sapply(unique_t,
                   function(x){
                     if(x %in% colnames(fm$Kappa)){
                       return(rmultinom(1, 1, fm$Kappa[,as.character(x)]))
                     } else {
                       last_kappa <- fm$Kappa[,ncol(fm$Kappa)]
                       steps <- x - as.numeric(colnames(fm$Kappa)[ncol(fm$Kappa)])
                       if(steps < 0){
                         stop("Backcasting not supported.")
                       }
                       new_kappa <- last_kappa %*% .mpower(fm$TransitionKernel, steps)
                       return(rmultinom(1, 1, new_kappa))
                     }}))
  colnames(states) <- unique_t
  states_ind <- states[,match(monad[,tid], colnames(fm$Kappa))]
  if(parametric_mm){
    a <- .pi.hat(X_m, fm$MonadCoef)
    pind <- lapply(a,
                   function(x){
                     apply(x, 2,
                           function(alpha){
                             res <- rgamma(n_blk, alpha, 1)
                             res/sum(res)
                           }
                     )})
    p <- .e.pi(pind, states_ind)
  } else {
    p <- fm$MixedMembership
  }

  z <- .getZ(p[,s_ind])
  w <- .getZ(p[,r_ind])

  eta_dyad <- array(0.0, n_dyad) 
  for(i in 1:n_dyad){ 
    eta_dyad[i] <- z[,i] %*% fm$BlockModel %*% w[,i]
  }
  eta_dyad <- eta_dyad + c(X_d %*% (fm$DyadCoef))

  probs <- plogis(eta_dyad)
  res <- rbinom(length(probs), 1, probs)
  return(res)
}
