#' Simulate a complete sociomatrix from an \code{mmsbmB} object
#'
#' The function generates one sample network from the posterior predictive of the model represented by a fitted \code{mmsbmB} object.
#'
#' @param object An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}
#' @param nsim Number of networks to simulate
#' @param seed RNG seed.
#' @param new.data.dyad An optional \code{data.frame} object. If not \code{NULL}, use these 
#'                      dyadic predictor values instead of those used to fit the original model.
#' @param new.data.monad1 An optional \code{data.frame} object. See \code{new.data.dyad}. 
#' @param new.data.monad2 An optional \code{data.frame} object. See \code{new.data.dyad}. 
#' @param ... Currently ignored
#' @return List of length \code{nsim} of simulated networks. 
#'         If \code{new.data.dyad = NULL}, each element is a vector of length \code{nrow(object$dyadic.data)}. 
#'         Else, vector of length \code{nrow(new.data.dyad)}. If \code{seed} is not NULL, return object
#'         includes its value as attribute "seed".
#' @method simulate mmsbmB
#'
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#'

simulate.mmsbmB <- function(object,
                            nsim = 1,
                            seed = NULL,
                            new.data.dyad = NULL,
                            new.data.monad1  = NULL,
                            new.data.monad2  = NULL,
                            parametric_mm = FALSE,
                            ...)
{
  if(!is.null(seed)){ set.seed(seed) }
  
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
    tid <- "(tid)" # add this line
    #sid <- object$forms$senderID
    #rid <- object$forms$receiverID
    if(is.null(object$forms$timeID)){ tid <- "(tid)"} else{tid <- object$forms$timeID}
    dyad <- object$dyadic.data
  }
  #Family 1 monadic
  if(!is.null(new.data.monad1)){
    nid1 <- object$forms$nodeID1
    monad1 <- new.data.monad1
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      monad1[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    if(is.null(object$forms$nodeID1)){nid1 <- "(nid)"
     tid <- "(tid)" # add this line
     }else{nid1<-object$forms$nodeID1
      tid <- "(tid)" # add this line
      }#nid1 <- "(nid1)"
    #tid <- "(tid)"
    monad1 <- object$monadic.data[[1]]
  }
  
  
  #Family 2 monadic
  if(!is.null(new.data.monad2)){
    nid2 <- object$forms$nodeID2
    monad2 <- new.data.monad2
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      monad2[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    if(is.null(object$forms$nodeID2)){nid2 <- "(nid)"
     tid <- "(tid)" # add this line
     }else{nid2<-object$forms$nodeID2
      tid <- "(tid)" # add this line
      }
    monad2 <- object$monadic.data[[2]]
  }
  
  if(! (tid %in% c(names(dyad), names(monad1)) |tid %in% c(names(dyad), names(monad2)) )){
    stop("Dynamic model estimated, but no timeID provided in new data.")
  }
  
  X_d <- model.matrix(eval(object$forms$formula.dyad), dyad)
  
  #Family 1 X_m
  if(is.null(object$forms$formula.monad1)){
    X_m1 <- model.matrix(~ 1, data = monad1)
  } else {
    X_m1 <- model.matrix(eval(object$forms$formula.monad1), monad1)
  }
  #Family 2 X_m
  if(is.null(object$forms$formula.monad2)){
    X_m2 <- model.matrix(~ 1, data = monad2)
  } else {
    X_m2 <- model.matrix(eval(object$forms$formula.monad2), monad2)
  }
  
  if(length(object$DyadCoef)==0){
    object$DyadCoef <- as.vector(0)
  } else {
    object$DyadCoef <- c(0, object$DyadCoef)
  }
  
  
  s_ind <- match(paste(dyad[,sid],dyad[,tid],sep="@"),
                 paste(monad1[,nid1],monad1[,tid],sep="@"))
  r_ind <- match(paste(dyad[,rid],dyad[,tid],sep="@"),
                 paste(monad2[,nid2],monad2[,tid],sep="@"))
  n_blk1 <- object$n_blocks1
  n_blk2 <- object$n_blocks2
  n_dyad <- nrow(X_d)
  
  unique_t <- unique(monad1[,tid]) #only for monad1 currently
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
    states<-t(states) #add this line
    colnames(states) <- unique_t # THIS IS CAUSING ERROR
    states_ind <- states[,match(monad1[,tid], colnames(object$Kappa))]#only for monad1 currently
    if(parametric_mm){
      #Family 1
      a <- .pi.hat(X_m1, object$MonadCoef1)
      pind <- lapply(a,
                     function(x){
                       apply(x, 2,
                             function(alpha){
                               res <- rgamma(n_blk1, alpha, 1)
                               res/sum(res)
                             }
                       )})
      p1 <- .e.pi(pind, states_ind)
      #Family 2
      a <- .pi.hat(X_m2, object$MonadCoef2)
      pind <- lapply(a,
                     function(x){
                       apply(x, 2,
                             function(alpha){
                               res <- rgamma(n_blk2, alpha, 1)
                               res/sum(res)
                             }
                       )})
      p2 <- .e.pi(pind, states_ind)
    } else {
      p1 <- object$MixedMembership1
      p2 <- object$MixedMembership2
    }
    
    z <- getZ(p1[,s_ind])
    w <- getZ(p2[,r_ind])
    
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