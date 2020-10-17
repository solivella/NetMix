#' Parametric bootstrap for the dynMMSBM
#'
#' The function performs a parametric bootstrap procedure to obtain uncertainty estimates for dynMMSBM coefficients.
#'
#' @param fm An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}. 
#' @param iter number of bootstrap iterations
#' @param parallel logical; indicates whether the bootstrap iterations should be run in parallel
#'
#'     
#' @return List of length \code{iter}.  Each entry contains a list of monadic1, monadic2 and dyadic coefficient estimates from the parametric bootstrap.
#
#' 


boot.mmsbmB <- function(fm, iter, parallel=TRUE){ 
  BootResClass <- function(res1=NULL,res2=NULL, res3=NULL){
    me <- list(res1 = res1,
               res2 = res2)
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }
  ctrl <- eval(fm$call$mmsbm.control)
  ctrl$phi1_init_t <- fm$`MixedMembership 1`
  ctrl$phi2_init_t <- fm$`MixedMembership 2`
  ##deal with if as.factor() name doesn't allow easy model.matrix formation with the formulas called from fm
  if(any(str_detect(names(fm$monadic1.data),pattern="as.factor"))){
    tmp<-as.data.frame(fm$monadic1.data[,which(str_detect(names(fm$monadic1.data),pattern="as.factor"))])
    tmp_name<-names(fm$monadic1.data)[str_detect(names(fm$monadic1.data),pattern="as.factor")]
    tmp_name<-gsub("as.factor\\(", "", tmp_name)
    tmp_name<-gsub("\\)", "", tmp_name)
    names(tmp)<-tmp_name
    fm$monadic1.data<-cbind(fm$monadic1.data,tmp)
  }
  if(any(str_detect(names(fm$monadic2.data),pattern="as.factor"))){
    tmp<-as.data.frame(fm$monadic2.data[,which(str_detect(names(fm$monadic2.data),pattern="as.factor"))])
    tmp_name<-names(fm$monadic2.data)[str_detect(names(fm$monadic2.data),pattern="as.factor")]
    tmp_name<-gsub("as.factor\\(", "", tmp_name)
    tmp_name<-gsub("\\)", "", tmp_name)
    names(tmp)<-tmp_name
    fm$monadic2.data<-cbind(fm$monadic2.data,tmp)
  }
  
  if(parallel){
    require(foreach)
    require(doMC)
    registerDoMC(parallel::detectCores()/2)
    boot.res <- foreach(1:iter) %dopar% {
      fm$dyadic.data$sim <- predict.mmsbmB(fm, type="response", posterior.pi=T) 
      fit <- mmsbmB(formula.dyad = update.formula(fm$call[[2]], sim ~ .), 
                   formula.monad1 = as.formula(fm$call[[3]]),
                   formula.monad2 = as.formula(fm$call[[4]]),
                   senderID = "(sid)", receiverID = "(rid)",
                   timeID = "(tid)", nodeID1 = "(nid1)", nodeID2 = "(nid2)",
                   data.dyad=fm$dyadic.data, data.monad1=fm$monadic1.data,data.monad2=fm$monadic2.data,
                   n.blocks1 = fm$call$n.blocks1, n.blocks2 = fm$call$n.blocks2,
                   n.hmmstates = fm$call$n.hmmstates,nodes2 = fm$forms$nodes2,
                   npred2 = fm$forms$npred2, 
                   directed=fm$call$directed,
                   mmsbm.control = ctrl)
      perm.vecs1 <- clue::solve_LSAP(log(fit$`MixedMembership 1` ) %*% log(t(fm$`MixedMembership 1` )), maximum=T)
      perm.vecs2 <- clue::solve_LSAP(log(fit$`MixedMembership 2` ) %*% log(t(fm$`MixedMembership 2` )), maximum=T)
      kap.vecs <- clue::solve_LSAP(log(fit$Kappa) %*% log(t(fm$Kappa)), maximum=T)
      result <- BootResClass()
      result$res1 <- fit$MonadCoef1[,perm.vecs1,kap.vecs]
      result$res2 <- fit$MonadCoef2[,perm.vecs2,kap.vecs]
      result$res3 <- fit$DyadCoef
      return(result)
    }
  }
  
  for(i in 1:iter){
    fm$dyadic.data$sim <- predict(fm, type="response", posterior.pi=T) 
    fit <- mmsbmB(formula.dyad = update.formula(fm$call[[2]], sim ~ .), 
                  formula.monad1 = as.formula(fm$call[[3]]),
                  formula.monad2 = as.formula(fm$call[[4]]),
                  senderID = "(sid)", receiverID = "(rid)",
                  timeID = "(tid)", nodeID1 = "(nid1)", nodeID2 = "(nid2)",
                  data.dyad=fm$dyadic.data, data.monad1=fm$monadic1.data,data.monad2=fm$monadic2.data,
                  n.blocks1 = fm$call$n.blocks1, n.blocks2 = fm$call$n.blocks2,
                  n.hmmstates = fm$call$n.hmmstates,nodes2 = fm$forms$nodes2,
                  npred2 = fm$forms$npred2, 
                  directed=fm$call$directed,
                  mmsbm.control = ctrl)
    perm.vecs1 <- clue::solve_LSAP(log(fit$`MixedMembership 1` ) %*% log(t(fm$`MixedMembership 1` )), maximum=T)
    perm.vecs2 <- clue::solve_LSAP(log(fit$`MixedMembership 2` ) %*% log(t(fm$`MixedMembership 2` )), maximum=T)
    kap.vecs <- clue::solve_LSAP(log(fit$Kappa) %*% log(t(fm$Kappa)), maximum=T)
    
    result <- BootResClass()
    result$res1 <- fit$MonadCoef1[,perm.vecs1,kap.vecs]
    result$res2 <- fit$MonadCoef2[,perm.vecs2,kap.vecs]
    result$res3 <- fit$DyadCoef
    return(result)
  }
}

