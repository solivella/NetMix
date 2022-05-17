#' Predict edges based on estimated mmsbmB model
#'
#' The function produces averaged predicted mixed membership plot for a given monadic covariate 
#' Defaults into using (1) min-to-max range of continuous covariate, (2) all levels of factor covariate
#'
#' @param object Object of class \code{mmsbmB}.
#' @param cov \code{string} describing covariate, found in data.monad.
#' @param data.monad \code{data.frame} object. 
#' @param family value taking 1 or 2 for family.   
#'     
#' @return \code{ggplot} object that can be further modified: y axis is average predicted MM, x axis covariate values, plot broken down by latent groups.
#'  
#' @author Kosuke Imai (imai@@harvard.edu), Adeline Lo (aylo@@wisc.edu), Santiago Olivella (olivella@@unc.edu), Tyler Pratt (tyler.pratt@@yale.edu) 
#' 
#' @example tests/Examples/cosponsorship.rmd
#' 
covMM<-function(object,cov,data.monad,family=1,...)
{
  require(stringr)
  require(ggplot2)
  blk<-ifelse(family==1,object$n_blocks1,object$n_blocks2)
  #checks
  if(!cov%in%colnames(data.monad)){ stop("cov must be a variable in data.monad")}
  if(class(data.monad[,cov])!="numeric"&class(data.monad[,cov])!="integer"&class(data.monad[,cov])!="factor"){
    stop("cov must be a factor, numeric or integer variable in data.monad")
  }
  if(class(data.monad[,cov])=="numeric"|class(data.monad[,cov])=="integer"){
    vals<-seq(min(data.monad[,cov]),max(data.monad[,cov]),length.out=nrow(data.monad))
    tmp<-lapply(vals,function(x){
      data.monad[,cov]<-x
      if(family==1){predict.mmsbmB(object,new.data.monad1=data.monad,type="mm")
      }else{ predict.mmsbmB(object,new.data.monad2=data.monad,type="mm")}
    })
    #get first element in list (mm for family 1/2)
    if(family==1){tmp2<-lapply(tmp,"[[",1)}else{tmp2<-lapply(tmp,"[[",2)}
    #take the mean of elements of list of matrices
    tmp3<-do.call("cbind",lapply(tmp2,function(x){rowMeans(x)}))
    #create plot data
    x <- rep(vals,times=blk)# x Axis: ideolvalue
    value <- as.vector(t(tmp3)) # y Axis: MM by group (unroll/row)
    group <- rep(paste("Group",1:blk,sep=" "),each=length(vals))# group, one shape per group
    data <- data.frame(x, value, group)
    p<-ggplot(data, aes(x=x, y=value, fill=group)) + ylab("Avg Predicted Mixed Membership") + geom_area(alpha=0.6 , size=1, colour="black") 
  }else{#factor
    vals<-unique(data.monad[,cov])
    tmp<-lapply(vals,function(x){
      data.monad[,cov]<-factor(x,levels=levels(data.monad[,cov]))
      if(family==1){predict.mmsbmB(object,new.data.monad1=data.monad,type="mm")
      }else{ predict.mmsbmB(object,new.data.monad2=data.monad,type="mm")}})
    #get first element in list (mm for family 1/2)
    if(family==1){tmp2<-lapply(tmp,"[[",1)}else{tmp2<-lapply(tmp,"[[",2)}
    #take the mean of elements of list of matrices
    tmp3<-do.call("cbind",lapply(tmp2,function(x){rowMeans(x)}))
    #create plot data
    x <- rep(vals,times=blk)# x Axis
    value <- as.vector(t(tmp3)) # y Axis: MM by group (unroll/row)
    group <- rep(paste("Group",1:blk,sep=" "),each=length(vals)) # group, one shape per group
    data <- data.frame(x, value, group)
    p<-ggplot(data, aes(x=x, y=value, fill=group)) +  ylab("Avg Predicted Mixed Membership") + geom_bar(alpha=0.6 ,position="stack", stat="identity", color="black")
  }
  return(p)
}#end fn covMM
