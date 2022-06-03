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
#' @param samples Either NULL (default) or an integer number of predictions based on
#' \code{samples} vectors of regression parameters obtained from their asymptotic sampling distributions. 
#' @param type Character string. The default is to use the linear predictor of edges. The alternative
#'     "response" returns predicted edge probabilities, and "mm" returns predicted mixed-memberships.    
#'     
#' @return Vector of predicted edge probabilities (if \code{type="response"}), vector of linear
#' edge predictors (if \code{type="link"}), or matrix of predicted mixed-membership probabilities (otherwise). If \code{samples!=NULL},
#' one more dimension is added the return object, with one prediction per sample; for example, if \code{type="response"} and
#' \code{samples=100}, the return object is a matrix with 100 columns.  
#'  
#' @author Kosuke Imai (imai@@harvard.edu), Adeline Lo (aylo@@wisc.edu), Santiago Olivella (olivella@@unc.edu), Tyler Pratt (tyler.pratt@@yale.edu) 
#' 
#' @example tests/Examples/cosponsorship.rmd
#' 
predict.mmsbmB <- function(object,
                           new.data.dyad = NULL,
                           new.data.monad1  = NULL, 
                           new.data.monad2  = NULL, 
                           forecast = FALSE,
                           samples = 1,
                           type = c("link", "response", "mm"),
                           ...)
{
  require(stringr)
  type <- match.arg(type)
  if(samples>1 & (samples%%1 > 0)){
      stop("If not NULL, samples must be a positive integer.")
  }
  if(samples>1 & (forecast == TRUE)){
    stop("Multiple samples for forecasting models not implemented yet.")
  }
  #Set up new dyadic data
  if(!is.null(new.data.dyad)){
    sid <- object$forms$senderID
    rid <- object$forms$receiverID
    dyad <- new.data.dyad
    if(is.null(object$forms$timeID)){
      tid <- object$forms$timeID
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
  dform <- object$forms$formula.dyad
  if(any(grepl("missing", names(object$DyadCoef)))){
    dform <- update(as.formula(dform), 
                    paste(c("~ .", names(object$DyadCoef)[grep("missing", 
                                                               names(object$DyadCoef))]), collapse=" + "))
  }
  X_d <- model.matrix(eval(dform), dyad)
  if(length(object$DyadCoef)==0|length(labels(terms(object$forms$formula.dyad)))==0){
    object$DyadCoef <- as.vector(0)
  } else {
    object$DyadCoef <- c(0, object$DyadCoef)
  }
  #Monadic
  #Family 1
  if(!is.null(new.data.monad1)){
    nid1 <- ifelse(object$forms$nodeID[[1]] %in% colnames(new.data.monad1), object$forms$nodeID[[1]], "(nid)")
    monad1 <- new.data.monad1
    if(is.null(object$forms$timeID)){
      tid1 <- "(tid)"
      monad1[,tid1] <- 1
    } else {
      tid1 <- object$forms$timeID
    }
    C_mat1 <- matrix(0, ncol = object$forms$n.blocks[1], nrow = nrow(monad1))
  } else {
    nid1 <- "(nid)"
    tid1 <- "(tid)"
    monad1 <- object$monadic.data[[1]]
    C_mat1 <- object$CountMatrix1
  }
  if(any(str_detect(names(monad1),pattern="as.factor"))){
    #add vars names that are as.factor, without as.factor heading
    tmp<-as.data.frame(monad1[,which(str_detect(names(monad1),pattern="as.factor"))])
    tmp_name<-names(monad1)[str_detect(names(monad1),pattern="as.factor")]
    tmp_name<-gsub("as.factor\\(", "", tmp_name)
    tmp_name<-gsub("\\)", "", tmp_name)
    names(tmp)<-tmp_name
    monad1<-cbind(monad1,tmp)
  }
  #Family 2
  if(!is.null(new.data.monad2)){
    nid2 <- ifelse(object$forms$nodeID[[2]] %in% colnames(new.data.monad2), object$forms$nodeID[[2]], "(nid)")
    monad2 <- new.data.monad2
    if(is.null(object$forms$timeID)){
      tid2 <- "(tid)"
      monad2[,tid2] <- 1
    } else {
      tid2 <- object$forms$timeID
    }
    C_mat2 <- matrix(0, ncol = object$forms$n.blocks[2], nrow = nrow(monad2))
  } else {
    nid2 <- "(nid)"
    tid2 <- "(tid)"
    monad2 <- object$monadic.data[[2]]
    C_mat2 <- object$CountMatrix2
  }
  if(any(str_detect(names(monad2),pattern="as.factor"))){
    #add vars names that are as.factor, without as.factor heading
    tmp <- as.data.frame(monad2[,which(str_detect(names(monad2),pattern="as.factor"))])
    tmp_name <- names(monad2)[str_detect(names(monad2),pattern="as.factor")]
    tmp_name <- gsub("as.factor\\(", "", tmp_name)
    tmp_name <- gsub("\\)", "", tmp_name)
    names(tmp) <- tmp_name
    monad2 <- cbind(monad2,tmp)
  }
  
  #Blks, Formulas
  n_blk1 <- object$n_blocks1
  n_blk2 <- object$n_blocks2
  mform1 <- object$forms$formula.monad[[1]]
  mform2 <- object$forms$formula.monad[[2]]
  if(any(grepl("missing", rownames(object$MonadCoef1)))){
    mform1 <- update(as.formula(mform1), 
                     paste(c("~ .", rownames(object$MonadCoef1)[grep("missing", rownames(object$MonadCoef1))]), collapse=" + "))
  }
  if(any(grepl("missing", rownames(object$MonadCoef2)))){
    mform2 <- update(as.formula(mform2), 
                     paste(c("~ .", rownames(object$MonadCoef2)[grep("missing", rownames(object$MonadCoef2))]), collapse=" + "))
  }
  if(is.null(mform1)){
    X1_m <- model.matrix(~ 1, data = monad1)
  } else {
    X1_m <- model.matrix(eval(mform1), monad1)
  }
  if(is.null(mform2)){
    X2_m <- model.matrix(~ 1, data = monad2)
  } else {
    X2_m <- model.matrix(eval(mform2), monad2)
  }
  if(samples==1){
    alpha1 <- list(.compute.alpha(X1_m, object$MonadCoef1))
    alpha2 <- list(.compute.alpha(X2_m, object$MonadCoef2))
  } else {
      MonadCoef1_samples  <- MASS::mvrnorm(samples,
                                           c(object$MonadCoef1),
                                           object$vcov_monad1)
      alpha1 <- apply(MonadCoef1_samples, 1,
                      function(betas){
                          betas <- array(betas, dim(object$MonadCoef1),
                                         dimnames = dimnames(object$MonadCoef1))
                          .compute.alpha(X1_m, betas)
                      })
      MonadCoef2_samples  <- MASS::mvrnorm(samples,
                                           c(object$MonadCoef2),
                                           object$vcov_monad2)
      alpha2 <- apply(MonadCoef2_samples, 1,
                      function(betas){
                          betas <- array(betas, dim(object$MonadCoef2),
                                         dimnames = dimnames(object$MonadCoef2))
                          .compute.alpha(X2_m, betas)
                      })
     
  }
  #Produce p1, p2
  
  if(forecast){
    ts1 <- unique(monad1[,tid])
    ts2 <- unique(monad2[,tid])
    new_kappa <- as.matrix(object$Kappa[,ncol(object$Kappa)] %*% .mpower(object$TransitionKernel, forecast))
    new_kappa1 <- matrix(new_kappa, nrow=ncol(new_kappa), ncol=nrow(monad1[monad1[,tid]==ts1[1],]),byrow=FALSE)
    new_kappa2 <- matrix(new_kappa, nrow=ncol(new_kappa), ncol=nrow(monad2[monad2[,tid]==ts2[1],]),byrow=FALSE)
    if(length(ts) > 1){
      for(t in 2:length(ts)){
        new_kappa <- rbind(new_kappa, new_kappa[t-1,] %*% .mpower(object$TransitionKernel, forecast))
        new_kappa1 <- cbind(new_kappa1, matrix(new_kappa[t,], nrow=ncol(new_kappa), ncol=nrow(monad1[monad1[,tid]==ts1[t],]),byrow=FALSE))
        new_kappa2 <- cbind(new_kappa2, matrix(new_kappa[t,], nrow=ncol(new_kappa), ncol=nrow(monad2[monad2[,tid]==ts2[t],]),byrow=FALSE))
      }
    }
    p1 <- .e.pi(alpha1, new_kappa1, C_mat1)
    p2 <- .e.pi(alpha2, new_kappa2, C_mat2)
  } else {
    #if(!(tid %in% colnames(monad1))){tid <- "(tid)"}
    p1 <- vapply(seq.int(length(alpha1)),
           function(x){
             .e.pi(alpha1[[1]], object$Kappa[,as.character(monad1[,tid])], C_mat1)
           },
           array(0, dim(alpha1[[1]][[1]]), dimnames = dimnames(alpha1[[1]][[1]])))
    p2 <- vapply(seq.int(length(alpha2)),
                 function(x){
                   .e.pi(alpha2[[1]], object$Kappa[,as.character(monad2[,tid])], C_mat2)
                 },
                 array(0, dim(alpha2[[1]][[1]]), dimnames = dimnames(alpha2[[1]][[1]])))
  }
  
  if(type=="mm"){
    return(list(MixedMembership1=p1, MixedMembership2=p2))
  }
  
  #Produce pi
  if(!sid%in%colnames(dyad)){tmp_sid <- object$forms$senderID}else{tmp_sid<-sid}
  if(!rid%in%colnames(dyad)){tmp_rid <- object$forms$receiverID}else{tmp_rid<-rid}
  if(!tid%in%colnames(dyad)){tmp_tid <- object$forms$timeID}else{tmp_tid<-tid}
  
  if(tmp_tid=="(tid)"&tid1!="(tid)"&any(names(monad1)%in%"(tid)")){
    temp_tid<-"(tid)"
    s_ind <- match(paste(dyad[,tmp_sid],dyad[,tmp_tid],sep="@"), 
                   paste(monad1[,nid1],monad1[,temp_tid],sep="@"))
  }else{
    s_ind <- match(paste(dyad[,tmp_sid],dyad[,tmp_tid],sep="@"), 
                   paste(monad1[,nid1],monad1[,tid1],sep="@"))
  }
  
  if (tmp_tid=="(tid)"&tid2!="(tid)"&any(names(monad2)%in%"(tid)")){
    temp_tid<-"(tid)"
    r_ind <- match(paste(dyad[,tmp_rid],dyad[,tmp_tid],sep="@"), 
                   paste(monad2[,nid2],monad2[,temp_tid],sep="@"))
  }else{
  r_ind <- match(paste(dyad[,tmp_rid],dyad[,tmp_tid],sep="@"), 
                 paste(monad2[,nid2],monad2[,tid2],sep="@"))
  }
  pi_s <- p1[,s_ind, ,drop=FALSE]
  pi_r <- p2[,r_ind, ,drop=FALSE]
  
  n_dyad <- nrow(X_d)
  eta_dyad <- array(0.0, c(n_dyad, ifelse(samples==1,1,samples)))
  for(j in 1:ncol(eta_dyad)){
    for(i in 1:n_dyad){ 
      eta_dyad[i,j] <- pi_s[,i,j] %*% object$BlockModel %*% pi_r[,i,j] #predicted
    }
  }
  if(samples>1){
    dyad_gamma <-  t(MASS::mvrnorm(samples,
                                 object$DyadCoef[-1],
                                 object$vcov_dyad))
  } else {
    dyad_gamma <-  c(object$DyadCoef[-1])
  }
  eta_dyad <- eta_dyad + (X_d[,-1] %*% dyad_gamma)#add dyadic
  if(type=="link"){
    res<-eta_dyad
  } else {
    res<-plogis(eta_dyad)#return via logit link
  }#evaluate prediction memberships of senators/bills; return pis and find max group for each senator/bill under predicted scenario
  return(res)
}
