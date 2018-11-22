simulate.mmsbm <- function(fm, 
                           new.data.dyad = NULL,
                           new.data.monad  = NULL, 
                           parametric_mm = FALSE)
{
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
  if(!is.null(new.data.monad)){
    nid <- fm$call$nodeID
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
  if(! (tid %in% c(names(dyad), names(monad)))){
    stop("Dynamic model estimated, but no timeID provided in new data.")
  }
  
  X_d <- model.matrix(eval(fm$call$formula.dyad), dyad)
  if(is.null(fm$call$formula.monad)){
    X_m <- model.matrix(~ 1, data = monad)
  } else {
    X_m <- model.matrix(eval(fm$call$formula.monad), monad)
  }
  if(length(fm$DyadCoef)==0){
    fm$DyadCoef <- as.vector(0)
  }
  s_ind <- match(dyad[,sid], monad[,nid])
  r_ind <- match(dyad[,rid], monad[,nid])
  t_ind <- match(dyad[,tid], monad[,tid])
  n_blk <- fm$n_blocks
  
  unique_t <- unique(monad[,tid])
  states <- as.matrix(sapply(unique_t,
                   function(x){
                     if(x %in% colnames(fm$Kappa)){
                       return(rmultinom(1, 1, fm$Kappa[,x]))
                     } else {
                       last_kappa <- fm$Kappa[,ncol(fm$Kappa)]
                       steps <- x - as.numeric(colnames(fm$Kappa)[ncol(fm$Kappa)])
                       if(steps < 0){
                         stop("Backcasting not supported.")
                       }
                       new_kappa <- last_kappa %*% .mpower(fm$TransitionKernel, steps)
                       return(rmultinom(1, 1, new_kappa))
                     }}))
  names(states) <- unique_t
  states_ind <- states[,match(monad[,tid], colnames(fm$Kappa))]
  if(parametric_mm){
    ph <- .pi.hat(X_m, fm$MonadCoef)
    a <- lapply(ph,
                function(x){
                  x * fm$MMConcentration
                })
    pind <- lapply(a,
                   function(x){
                     apply(x, 2,
                           function(alpha){
                             prop.table(rgamma(n_blk, alpha, 1))
                           }
                     )})
    p <- .e.pi(pind, states_ind)
  } else {
    p <- fm$MixedMembership
  }
  z <- getZ(p[,s_ind])
  w <- getZ(p[,r_ind])
  if(length(fm$DyadCoef) >= 1){
    eta_dyad <- X_d[,-1] %*% t(fm$DyadCoef)
  } else {
    eta_dyad <- X_d %*% t(fm$DyadCoef)
  }
  for(a in 1:n_blk){ 
    for(b in 1:n_blk){
      eta_dyad <- eta_dyad + (z[a,]*w[b,]*fm$BlockModel[a,b])
    }
  }
  probs <- plogis(eta_dyad)
  return(rbinom(length(probs), 1, probs))
}