predict.mmsbm <- function(fm, 
                          new.data.dyad = NULL,
                          new.data.monad  = NULL, 
                          type = c("prediction","expectation"), 
                          parametric_mm = FALSE,
                          out.sample = FALSE,
                          outcome = c("probability", "response")){
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
  X_d <- model.matrix(eval(fm$call$formula.dyad), dyad)[,-1]
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
  if(!out.sample){
    if(!parametric_mm){
      p <- fm$MixedMembership
    } else {
      ph <- .pi.hat(X_m, fm$MonadCoef)
    }
    eta_dyad <- X_d %*% as.matrix(fm$DyadCoef)
    if(type=="expectation"){
      if(parametric_mm){
        p <- .e.pi(ph, fm$Kappa[,t_ind])
      }
      for(a in 1:n_blk){ 
        for(b in 1:n_blk){
          eta_dyad <- eta_dyad + (p[a,s_ind]*p[b,r_ind]*fm$BlockModel[a,b])
        }
      }
    } else { #prediction
      if(parametric_mm){
        s <- apply(fm$Kappa, 2, function(x){rmultinom(1, 1, x)})
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
        p <- .e.pi(pind, s)
      }
      z <- apply(p[,s_ind], 2, function(x){rmultinom(1, 1, x)})
      w <- apply(p[,r_ind], 2, function(x){rmultinom(1, 1, x)}) 
      for(a in 1:n_blk){ 
        for(b in 1:n_blk){
          eta_dyad <- eta_dyad + (z[a,]*w[b,]*fm$BlockModel[a,b])
        }
      }
    }
    probs <- plogis(eta_dyad)
    if(outcome=="probability"){
      return(probs)
    } else {
      return(rbinom(length(probs), 1, probs))
    }
    ## END OF REVISED VERSION
  } else {
    tid <- ifelse(all(colnames(fm$monadic.data) %in% colnames(monad)), "(tid)", fm$call$timeID)
    kappa.last <- fm$Kappa[,ncol(fm$Kappa)]
    kappa.new <- list()
    kappa.new[[1]] <- fm$TransitionKernel %*% kappa.last
    if(length(unique(dyad$year))>1){
      for(i in 2:length(unique(dyad$year))){
        kappa.new[[i]] <- fm$TransitionKernel %*% kappa.new[[i-1]]
      }
    }
    new.pis <- list()
    for(i in 1:fm$n_states){
      ifelse(!is.null(fm$call$formula.monad),
             monad.vars <- cbind(rep(1, nrow(monad)), get_all_vars(fm$call$formula.monad, monad)),
             monad.vars <- rep(1, nrow(fm$monadic.data)))
      new.pis[[i]] <- as.matrix(monad.vars) %*% fm$MonadCoef[,,i]
    }
    mean.pis <- new.pis[[1]]
    for(j in 1:length(unique(monad[,tid]))){
      kappa.rel <- kappa.new[[j]]
      mean.pi.sub <- mean.pis[monad[,tid]==unique(monad[,tid])[j],]
      for(q in 1:nrow(mean.pi.sub)){
        pi.mean <- rep(0, ncol(mean.pi.sub))
        for(v in 1:fm$n_states){
          pi.mean <- pi.mean + new.pis[[v]][q,] * kappa.new[[j]][v]
        }
        mean.pi.sub[q,] <- pi.mean
      }
      mean.pis[monad[,tid]==unique(monad[,tid])[j],] <- mean.pi.sub
    }
    pi1 <- t(mean.pis[match(dyad[,fm$call$senderID], monad[,fm$call$nodeID]),])
    pi2 <- t(mean.pis[match(dyad[,fm$call$receiverID], monad[,fm$call$nodeID]),])
  }
  
  z1 <- apply(pis[,match(dyad[,sid], monad[,nid])], 2, function(w){rmultinom(1, 1, w)})
  z2 <- apply(pis[,match(dyad[,rid], monad[,nid])], 2, function(w){rmultinom(1, 1, w)})
  
  Z <- model.matrix(as.formula(fm$call$formula.dyad), data = dyad)
  est.ties <- Z %*% t(fm$DyadCoef)
  
  for(a in 1:nrow(pi1)){ 
    for(b in 1:nrow(pi2)){
      est.ties <- est.ties + (pi1[a,]*pi2[b,]*fm$BlockModel[a,b])
    }
  }
  if(type=="probability"){
    return(exp(est.ties) / (1 + exp(est.ties)))}
  if(type=="response"){
    return(rbinom(length(est.ties), 1, prob=exp(est.ties) / (1 + exp(est.ties))))
  }
}

