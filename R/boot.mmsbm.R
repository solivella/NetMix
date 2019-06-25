#' Parametric bootstrap for the dynamic mixed-membership SBM
#'
#' The function performs a parametric bootstrap procedure on \code{mmsbm} models.
#'
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}. 
#' @param iter number of bootstrap iterations
#' @param level signficance level for bootstrap confidence intervals
#' @param full_obj logical; return full \code{mmsbm} objects, or bootstrapped CI for regression coefficients?
#' @param parallel logical; indicates whether the bootstrap iterations should be run in parallel
#'
#'     
#' @return If \code{full_obj=TRUE}, list of length \code{iter}. Each entry contains an object of class \code{mmsbm}. If \code{full_obj=FALSE}, a named list:
#' \describe{
#'       \item{MonadCoef}{List with as many elements as there are HMM states in \code{fm}. Each list is an array with as many rows as there are monadic predictors,
#'       and two columns (an upper and lower bound for the corresponding CI)}
#'       \item{DyadCoef}{Array with as many rows as there are dyadic predictors,
#'       and two columns (an upper and lower bound for the corresponding CI)}
#'     }
#' 


boot.mmsbm <- function(fm, iter, level = 0.9, full_obj = FALSE, parallel = TRUE){ 
  BootResClass <- function(monadic=NULL,dyadic=NULL){
    me <- list(monadic = monadic,
               dyadic = dyadic)
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }
  
  ctrl_tmp <- fm$forms$mmsbm.control
  ctrl_tmp$threads <- 1
  ctrl_tmp$verbose <- FALSE
  
  result <- vector("list", iter) 
  require(foreach)
  `%our_do%` <- if (parallel) `%dopar%` else `%do%`
  if(parallel){
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
                 data.dyad = new_dyad,
                 data.monad = monad.data,
                 n.blocks = nb,
                 n.hmmstates = ns,
                 directed = direct,
                 mmsbm.control = ctrl_tmp)
    loss.mat.phi<- fit$MixedMembership %*% t(fm$MixedMembership)
    phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
    orig_mm_names <- rownames(fit$MixedMembership)
    fit$MixedMembership <- fit$MixedMembership[phi_ord, ]
    rownames(fit$MixedMembership) <- orig_mm_names
    orig_g_names <- dimnames(fit$BlockModel)
    fit$BlockModel <- fit$BlockModel[phi_ord, phi_ord]
    dimnames(fit$BlockModel) <- orig_g_names
    if(fm$n_states > 1){
      loss.mat.kappa <- fit$Kappa %*% t(fm$Kappa)
      kappa_ord <- clue::solve_LSAP(t(loss.mat.kappa), TRUE)
    } else {
      kappa_ord <- 1
    }
    orig_coef_names <- dimnames(fit$MonadCoef)
    fit$MonadCoef <- fit$MonadCoef[, phi_ord, kappa_ord, drop = FALSE]
    dimnames(fit$MonadCoef) <- orig_coef_names
    orig_k_names <- rownames(fit$Kappa)
    fit$Kappa <- fit$Kappa[kappa_ord, , drop = FALSE]
    rownames(fit$Kappa) <- orig_k_names
    orig_tk_names <- dimnames(fit$TransitionKernel)
    fit$TransitionKernel <- fit$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
    dimnames(fit$TransitionKernel) <- orig_tk_names
    result <- BootResClass(fit$MonadCoef, fit$DyadCoef)#
    if(full_obj){
      return(fit)
    } else {
      return(result)
    }
  }
  if(!full_obj){
    monadic_coef <- vector("list", fm$n_states)
    for(i in 1:fm$n_states){
      all_coef_monad <- sapply(boot.res,
                               function(x){
                                 res <- c(NA, x$monadic[,,i, drop = FALSE])
                                 names(res) <- c(NA,apply(do.call(expand.grid, dimnames(x$monadic[,,i, drop = FALSE])[c(1,2)]), 1, paste, collapse="_"))
                                 return(res)
                               })
      monadic_coef[[i]] <- t(apply(all_coef_monad, 1, quantile, probs = c((1-level)/2, 0.5, 1 - (1-level)/2), na.rm=TRUE))
      monadic_coef[[i]] <- as.array(monadic_coef[[i]][complete.cases(monadic_coef[[i]]), ,drop = FALSE])
    }
    names(monadic_coef) <- paste("State",1:fm$n_states)
    all_coef_dyad <- sapply(boot.res,
                            function(x){
                              c(NA,x$dyadic)
                            })
    all_coef_dyad <- t(apply(all_coef_dyad, 1, quantile, probs = c((1-level)/2, 0.5, 1 - (1-level)/2), na.rm=TRUE))
    all_coef_dyad <- as.array(all_coef_dyad[complete.cases(all_coef_dyad),, drop=FALSE])
    boot.res <- list(MonadCoef = monadic_coef,
                     DyadCoef = all_coef_dyad)
  }
  return(boot.res)
}



