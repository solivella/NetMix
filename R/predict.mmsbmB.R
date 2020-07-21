#' Predict edges based on estimated mmsbmB model
#'
#' The function produces expected posterior edges based  
#' on estimated parameters and (optionally new) predictor data 
#'
#' @param fm Object of class \code{mmsbmB}.
#' @param new.data.dyad An optional \code{data.frame} object. 
#' @param new.data.monad1 An optional \code{data.frame} object. 
#' @param new.data.monad2 An optional \code{data.frame} object. 
#' @param forecast Boolean. Defaults to FALSE.
#' @param type Character string. The default is to use the linear predictor of edges. The alternative
#'     "response" returns predicted probabilities.    
#'     
#' @return If \code{new.data.dyad = NULL}, vector of length \code{nrow(fm$dyadic.data)}. Else, vector of length \code{nrow(new.data.dyad)}.
#'  
#' @author Kosuke Imai (imai@@harvard.edu), Adeline Lo (aylo@@wisc.edu), Santiago Olivella (olivella@@unc.edu), Tyler Pratt (tyler.pratt@@yale.edu) 
#' 
#' @example tests/Examples/cosponsorship.rmd
#' 
predict.mmsbmB <- function(object,
                           new.data.dyad = NULL,
                           new.data.monad1  = NULL, 
                           new.data.monad2  = NULL, 
                           parametric_mm = FALSE,
                           forecast = FALSE,
                           type = c("link", "response"),
                           ...)
{
  type <- match.arg(type)
  if(((!is.null(new.data.monad1) | !is.null(new.data.monad2)) | forecast) & !parametric_mm){
    stop("Must use parametric mixed-memberships when forecasting or when using new monadic data.")
  }  
  #Set up new dyadic data
  if(!is.null(new.data.dyad)){
    sid <- object$forms$senderID
    rid <- object$forms$receiverID
    dyad <- new.data.dyad
    if(is.null(object$forms$timeID)){
      #tid <- "(tid)"
      tid <- object$forms$timeID
      dyad[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    sid <- "(sid)" #ifelse(!is.null(object$forms$nodeID1),object$forms$nodeID1,"(sid)") 
    rid <- "(rid)" #ifelse(!is.null(object$forms$nodeID2),object$forms$nodeID2,"(rid)") 
    if(!is.null(new.data.monad1)|!is.null(new.data.monad2)){
      tid <-  object$forms$timeID  
    }else {tid <- "(tid)"}
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
  
  #Monadic
  #Family 1
  if(!is.null(new.data.monad1)){
    nid1 <- ifelse(object$forms$nodeID1 %in% colnames(new.data.monad1), object$forms$nodeID1, "(nid1)")
    monad1 <- new.data.monad1
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      monad1[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    #nid1 <- "(nid1)"
    nid1 <- object$forms$nodeID1
    #tid <- "(tid)"
    tid <- object$forms$timeID
    monad1 <- object$monadic1.data
  }
  #Family 2
  if(!is.null(new.data.monad2)){
    nid2 <- ifelse(object$forms$nodeID2 %in% colnames(new.data.monad2), object$forms$nodeID2, "(nid2)")
    monad2 <- new.data.monad2
    if(is.null(object$forms$timeID)){
      tid <- "(tid)"
      monad2[,tid] <- 1
    } else {
      tid <- object$forms$timeID
    }
  } else {
    #nid2 <- "(nid2)"
    nid2 <- object$forms$nodeID2
    #tid <- "(tid)"
    tid <- object$forms$timeID
    monad2 <- object$monadic2.data
  }
  
  #Blks, Formulas
  n_blk1 <- object$n_blocks1
  n_blk2 <- object$n_blocks2
  mform1 <- object$forms$formula.monad1
  mform2 <- object$forms$formula.monad2
  if(any(grepl("missing", rownames(object$MonadCoef1)))){
    mform1 <- update(as.formula(mform1), 
                     paste(c("~ .", rownames(object$MonadCoef1)[grep("missing", rownames(object$MonadCoef1))]), collapse=" + "))
  }
  if(any(grepl("missing", rownames(object$MonadCoef2)))){
    mform2 <- update(as.formula(mform2), 
                     paste(c("~ .", rownames(object$MonadCoef2)[grep("missing", rownames(object$MonadCoef2))]), collapse=" + "))
  }
  
  
  #Produce p1, p2
  if(!parametric_mm){
    p1 <- object$`MixedMembership 1` #push exp value pi return
    p2 <- object$`MixedMembership 2`
  } else {
    if(is.null(mform1)){ X1_m <- model.matrix(~ 1, data = monad1) } else { X1_m <- model.matrix(eval(mform1), monad1)}
    if(is.null(mform2)){ X2_m <- model.matrix(~ 1, data = monad2) } else { X2_m <- model.matrix(eval(mform2), monad2)}
    alpha1 <- .pi.hat(X1_m, object$MonadCoef1)
    alpha2 <- .pi.hat(X2_m, object$MonadCoef2)
    if(forecast){
      ts <- unique(monad1[,tid])
      new_kappa <- as.matrix(object$Kappa[,ncol(object$Kappa)] %*% .mpower(object$TransitionKernel, forecast))
      new_kappa1 <- matrix(new_kappa, nrow=ncol(new_kappa), ncol=nrow(monad1[monad1[,tid]==ts[1],]),byrow=FALSE)
      new_kappa2 <- matrix(new_kappa, nrow=ncol(new_kappa), ncol=nrow(monad2[monad2[,tid]==ts[1],]),byrow=FALSE)
      if(length(ts) > 1){
        for(t in 2:length(ts)){
          new_kappa <- rbind(new_kappa, new_kappa[t-1,] %*% .mpower(object$TransitionKernel, forecast))
          new_kappa1 <- cbind(new_kappa1, matrix(new_kappa[t,], nrow=ncol(new_kappa), ncol=nrow(monad1[monad1[,tid]==ts[t],]),byrow=FALSE))
          new_kappa2 <- cbind(new_kappa1, matrix(new_kappa[t,], nrow=ncol(new_kappa), ncol=nrow(monad2[monad2[,tid]==ts[t],]),byrow=FALSE))
        }
      }
      p1 <- .e.pi(lapply(alpha1, function(x)prop.table(x, 2)),
                  new_kappa1)
      p2 <- .e.pi(lapply(alpha2, function(x)prop.table(x, 2)),
                  new_kappa2)
    } else {
      if(!(tid %in% colnames(monad1))){tid <- "(tid)"}
      p1 <- .e.pi(lapply(alpha1, function(x)prop.table(x, 2)), #computes expected mm based on covariates
                  object$Kappa)
      # object$Kappa[,as.character(monad1[,tid])]) #needs changing to allow for multiple Kappa
      p2 <- .e.pi(lapply(alpha2, function(x)prop.table(x, 2)), 
                  object$Kappa)
      #[,as.character(monad2[,tid])])
    }
  }
  
  #Produce pi
  s_ind <- match(paste(dyad[,sid],dyad[,tid],sep="@"), 
                 paste(monad1[,nid1],monad1[,tid],sep="@"))
  r_ind <- match(paste(dyad[,rid],dyad[,tid],sep="@"), 
                 paste(monad2[,nid2],monad2[,tid],sep="@"))
  pi_s <- p1[,s_ind]
  pi_r <- p2[,r_ind]
  
  n_dyad <- nrow(X_d)
  eta_dyad <- array(0.0, n_dyad)
  for(i in 1:n_dyad){ 
    eta_dyad[i] <- pi_s[,i] %*% object$BlockModel %*% pi_r[,i] #predicted
  }
  eta_dyad <- eta_dyad + c(X_d %*% (object$DyadCoef))#add dyadic
  if(type=="link"){
    res<-c(eta_dyad)
  } else {
    res<-plogis(c(eta_dyad))#return via logit link
  }#evaluate prediction memberships of senators/bills; return pis and find max group for each senator/bill under predicted scenario
  return(res)
}