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
                           forecast = FALSE,
                           type = c("link", "response"),
                           ...)
{
  require(stringr)
  type <- match.arg(type)
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
    tmp<-as.data.frame(monad2[,which(str_detect(names(monad2),pattern="as.factor"))])
    tmp_name<-names(monad2)[str_detect(names(monad2),pattern="as.factor")]
    tmp_name<-gsub("as.factor\\(", "", tmp_name)
    tmp_name<-gsub("\\)", "", tmp_name)
    names(tmp)<-tmp_name
    monad2<-cbind(monad2,tmp)
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
  alpha1 <- .compute.alpha(X1_m, object$MonadCoef1)
  alpha2 <- .compute.alpha(X2_m, object$MonadCoef2)
  
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
    p1 <- .e.pi(alpha1, object$Kappa[,as.character(monad1[,tid])], C_mat1)
    p2 <- .e.pi(alpha2, object$Kappa[,as.character(monad2[,tid])], C_mat2)
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