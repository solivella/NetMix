#' Simulate a complete sociomatrix from an \code{mmsbm} object
#'
#' The function generates one sample network from the posterior predictive of the model represented by a fitted \code{mmsbm} object.
#'
#' @param object An object of class \code{mmsbm}, a result of a call to \code{mmsbm}
#' @param nsim Number of networks to simulate
#' @param seed RNG seed.
#' @param new.data.dyad An optional \code{data.frame} object. If not \code{NULL}, use these 
#'                      dyadic predictor values instead of those used to fit the original model.
#' @param new.data.monad An optional \code{data.frame} object. See \code{new.data.dyad}. 
#' @param parametric_mm Boolean. Should the variational posterior be used for sampling the mixed-memberships (\code{FALSE}), 
#'                      or should the mixed-meberships be formed using the parameters in the monadic regression equation (\code{TRUE})?
#'                      Defaults to \code{FALSE}. If \code{is.null(new.data.monad)=FALSE}, setting this to \code{FALSE} will produce an error.  
#' @param ... Currently ignored
#' @return List of length \code{nsim} of simulated networks. 
#'         If \code{new.data.dyad = NULL}, each element is a vector of length \code{nrow(object$dyadic.data)}. 
#'         Else, vector of length \code{nrow(new.data.dyad)}. If \code{seed} is not NULL, return object
#'         includes its value as attribute "seed". .
#' @method simulate mmsbm
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
#' ## Simulate 10 new networks
#' lazega_sim <- simulate(lazega_mmsbm, nsim = 10, seed = 123)
#' 

simulate.mmsbm <- function(object, 
                           nsim = 1,
                           seed = NULL,
                           new.data.dyad = NULL,
                           new.data.monad  = NULL, 
                           parametric_mm = FALSE,
                           ...)
{
  if(!is.null(seed)){
    set.seed(seed)
  }
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
  if(!is.null(new.data.monad)){
    nid <- object$forms$nodeID
    monad <- new.data.monad
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      monad[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    nid <- "(nid)"
    tid <- "(tid)"
    monad <- object$monadic.data
  }
  if(! (tid %in% c(names(dyad), names(monad)))){
    stop("Dynamic model estimated, but no timeID provided in new data.")
  }
  X_d <- model.matrix(object$forms$formula.dyad, dyad)
  if(is.null(object$forms$formula.monad)){
    X_m <- model.matrix(~ 1, data = monad)
  } else {
    X_m <- model.matrix(object$forms$formula.monad, monad)
  }
  if(length(object$DyadCoef)==0){
    object$DyadCoef <- as.vector(0)
  } else {
    object$DyadCoef <- c(0, object$DyadCoef)
  }
  s_ind <- match(paste(dyad[,sid],dyad[,tid],sep="@"), 
                 paste(monad[,nid],monad[,tid],sep="@"))
  r_ind <- match(paste(dyad[,rid],dyad[,tid],sep="@"), 
                 paste(monad[,nid],monad[,tid],sep="@"))
  n_blk <- object$n_blocks
  n_dyad <- nrow(X_d)
  
  unique_t <- unique(monad[,tid])
  res <- lapply(seq_len(nsim), function(n){
    states <- as.matrix(sapply(unique_t,
                               function(x){
                                 if(x %in% colnames(object$Kappa)){
                                   return(rmultinom(1, 1, object$Kappa[,as.character(x)]))
                                 } else {
                                   last_kappa <- object$Kappa[,ncol(object$Kappa)]
                                   steps <- x - as.numeric(colnames(object$Kappa)[ncol(object$Kappa)])
                                   if(steps < 0){
                                     stop("Backcasting not supported.")
                                   }
                                   new_kappa <- last_kappa %*% .mpower(object$TransitionKernel, steps)
                                   return(rmultinom(1, 1, new_kappa))
                                 }}))
    colnames(states) <- unique_t
    states_ind <- states[,match(monad[,tid], colnames(object$Kappa))]
    if(parametric_mm){
      a <- .pi.hat(X_m, object$MonadCoef)
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
      p <- object$MixedMembership
    }
    
    z <- .getZ(p[,s_ind])
    w <- .getZ(p[,r_ind])
    
    eta_dyad <- array(0.0, n_dyad) 
    for(i in 1:n_dyad){ 
      eta_dyad[i] <- z[,i] %*% object$BlockModel %*% w[,i]
    }
    eta_dyad <- eta_dyad + c(X_d %*% (object$DyadCoef))
    
    probs <- plogis(eta_dyad)
    res_int  <- rbinom(length(probs), 1, probs)
    return(res_int)
  })
  if(!is.null(seed)){
    attr(res, "seed") <- seed
  }
  return(res)
}
