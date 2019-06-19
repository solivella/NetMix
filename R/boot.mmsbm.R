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


boot.mmsbm <- function(fm, iter, parallel_boot=TRUE){ 
  BootResClass <- function(monadic=NULL,dyadic=NULL){
    me <- list(monadic = monadic,
               dyadic = dyadic)
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }
  
  ctrl_tmp <- fm$forms$mmsbm.control
  ctrl_tmp$phi_init_t <- fm$MixedMembership
  ctrl_tmp$kappa_init_t <- fm$Kappa
  ctrl_tmp$beta_init <- fm$MonadCoef
  ctrl_tmp$gamma_init <- fm$DyadCoef
  ctrl_tmp$threads <- 1
  ctrl_tmp$verbose <- FALSE
  
  result <- vector("list", iter) 
  require(foreach)
  `%our_do%` <- if (parallel_boot) `%dopar%` else `%do%`
  if(parallel_boot){
    require(doParallel)
    registerDoParallel()
  }
  sims <- replicate(iter, simulate(fm))
  dyad.formula <- update.formula(fm$forms$formula.dyad, Y~.)
  monad.formula <- fm$forms$formula.monad
  dyad.data <- fm$dyadic.data
  monad.data <- fm$monadic.data
  nb <- fm$n_blocks
  ns <- fm$n_states
  direct <- fm$directed
  boot.res <- foreach(i = 1:iter, .packages="NetMix") %our_do% {
    new_dyad <- dyad.data
    new_dyad$Y <- sims[,i]
    fit <- mmsbm(formula.dyad = dyad.formula,
                 formula.monad = monad.formula,
                 senderID = "(sid)", receiverID = "(rid)",
                 timeID = "(tid)", nodeID = "(nid)",
                 data.dyad=new_dyad,
                 data.monad=monad.data,
                 n.blocks = nb,
                 n.hmmstates = ns,
                 directed=direct,
                 mmsbm.control = ctrl_tmp)
    result <- BootResClass()#
    result$monadic <- fit$MonadCoef
    result$dyadic <- fit$DyadCoef
    return(result)
  }
  return(boot.res)
}



