#' Parametric bootstrap for the dynMMSBM
#'
#' The function performs a parametric bootstrap procedure to obtain uncertainty estimates for dynMMSBM coefficients.
#'
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}. 
#' @param iter number of bootstrap iterations
#' @param parallel logical; indicates whether the bootstrap iterations should be run in parallel
#'
#'     
#' @return List of length \code{iter}.  Each entry contains a list of monadic and dyadic coefficient estimates from the parametric bootstrap.
#
#' 


boot.mmsbm <- function(fm, iter, parallel=TRUE){ 
  BootResClass <- function(monadic=NULL,dyadic=NULL){
    me <- list(monadic = monadic,
               dyadic = dyadic)
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }
  
  ctrl <- eval(fm$call$mmsbm.control)
  ctrl$phi_init_t <- fm$MixedMembership
  
  result <- vector("list", iter) 
  if(parallel){
    require(foreach)
    require(doMC)
    registerDoMC(parallel::detectCores()/2)
    boot.res <- foreach(1:iter) %dopar% {
      fm$dyadic.data$sim <- simulate(fm) 
      fit <- mmsbm(formula.dyad = update.formula(fm$call[[2]], sim ~ .), 
                   formula.monad = as.formula(fm$call[[3]]),
                   senderID = "(sid)", receiverID = "(rid)",
                   timeID = "(tid)", nodeID = "(nid)",
                   data.dyad=fm$dyadic.data, data.monad=fm$monadic.data,
                   n.blocks = fm$n_blocks, 
                   n.hmmstates = fm$n_states,
                   directed=fm$directed,
                   mmsbm.control = ctrl)
      perm.vecs <- clue::solve_LSAP(log(fit$MixedMembership) %*% log(t(fm$MixedMembership)), maximum=T)
      kap.vecs <- clue::solve_LSAP(log(fit$Kappa) %*% log(t(fm$Kappa)), maximum=T)
      result <- BootResClass()
      result[[i]]
      result$monadic <- fit$MonadCoef[,perm.vecs,kap.vecs]
      result$dyadic <- fit$DyadCoef
      return(result)
    }
    return(boot.res)
  }
  
  if(!parallel){
    for(i in 1:iter){
      fm$dyadic.data$sim <- simulate(fm) 
      fit <- mmsbm(formula.dyad = update.formula(fm$call[[2]], sim ~ .), 
                 formula.monad = as.formula(fm$call[[3]]),
                 senderID = "(sid)", receiverID = "(rid)",
                 timeID = "(tid)", nodeID = "(nid)",
                 data.dyad=fm$dyadic.data, data.monad=fm$monadic.data,
                 n.blocks = fm$n_blocks, 
                 n.hmmstates = fm$n_states,
                 directed=fm$directed,
                 mmsbm.control = ctrl)
      perm.vecs <- clue::solve_LSAP(log(fit$MixedMembership) %*% log(t(fm$MixedMembership)), maximum=T)
      kap.vecs <- clue::solve_LSAP(log(fit$Kappa) %*% log(t(fm$Kappa)), maximum=T)
      result <- BootResClass()
      result$monadic <- fit$MonadCoef[,perm.vecs,kap.vecs]
      result$dyadic <- fit$DyadCoef
      return(result)
    }
    return(result)
  }
}

