#' Generate plots for latent group membership by covariate from estimated mmsbmB model
#'
#' The function generates plots that illustrate covariate levels for each latent group.
#'
#' @param fm An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}. 
#' @param cov Character string identifying the monadic covariate to be shifted.  
#' @param family number 1 or 2 identifying the node family from which to find the monadic covariate/latent groups  
#' @param groupassign takes values 'expected' or 'max', for method of assigning groups 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 

covGroupB<-function(fm, cov=NULL, family=NULL, groupassign="expected"){
  if (!requireNamespace("sm", quietly = TRUE)) {
    stop("Package \"sm\" needed to produce requested plot. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package \"viridis\" needed to produce requested plot. Please install it.",
         call. = FALSE)
  }
  require(sm)
  require(viridis)
  #Checks
  if(is.null(cov)){cat("Error: no covariate selected.")}
  if(is.null(family)){cat("Error: no node family selected -- please choose 1 or 2.")}
  if(family==1){
    covariate<-fm$monadic.data[[1]][,cov]
  }else{
    covariate<-fm$monadic.data[[2]][,cov]
  }
  #if(class(covariate)!="numeric"|class(covariate)!="integer"|class(covariate)!="factor"){
    #scat("Error: covariate must be numeric/integer/factor class.")}
  
  if(family==1){tmp_blk<-fm$n_blocks1}else{tmp_blk<-fm$n_blocks2}
  if(family==1){tmp_membership<-fm$`MixedMembership1`}else{tmp_membership<-fm$`MixedMembership2`}
  ## Univariate Density Estimate Comparisons
  if(class(covariate)=="numeric"|class(covariate)=="integer"){
    tmp_col<-viridis(tmp_blk,option="E")
    groups<-apply(tmp_membership,2,which.max)
    groups.label <- factor(groups, levels= 1:tmp_blk,
                           labels = 1:tmp_blk) 
    sm.density.compare(covariate, groups, col=tmp_col, xlab=paste(cov), lwd=3)
    legend("topright", levels(as.factor(groups)), fill=tmp_col)
  }
  if(class(covariate)=="factor"& groupassign=="max"){
    groups<-apply(tmp_membership,2,which.max)
    group<-as.factor(rep(1:tmp_blk,length(levels(covariate))))
    proportion <-c(prop.table(table(groups,covariate),1))
    level<-rep(levels(covariate),each=length(unique(group)))#repeat by number of groups
    data <- data.frame(group,proportion,level)
    pos<-"fill"#"dodge"/"stack"/"fill"
    tmp_plot<-ggplot(data, aes(fill=group, y=level, x=proportion)) + ggtitle("Group proportions in each covariate level") +
      geom_bar(position=pos, stat="identity") + theme_bw() + scale_fill_viridis(discrete = T, option = "E", alpha=0.75)
    return(tmp_plot)
  }
  if(class(covariate)=="factor"& groupassign=="expected"){
    tmp_d<-vector("list",length(levels(covariate)))
    for(i in 1:length(levels(covariate))){
      if(length(tmp_membership[,which(covariate==levels(covariate)[i])])==tmp_blk){
        tmp_d[[i]]<-tmp_membership[,which(covariate==levels(covariate)[i])]/table(covariate)[i]
        #tmp_d[[i]]<-matrix(tmp_membership[,which(covariate==levels(covariate)[i])]/table(covariate)[i]  )
        #colnames(tmp_d[[i]])<-colnames(tmp_membership)[which(covariate==levels(covariate)[i])]
      }else{
      tmp_d[[i]]<-apply(tmp_membership[,which(covariate==levels(covariate)[i])],1,sum)/table(covariate)[i]
      }
    }
    
    group<-as.factor(rep(1:tmp_blk,length(levels(covariate))))
    proportion <-unlist(tmp_d)
    level<-rep(levels(covariate),each=length(unique(group)))#repeat by number of groups
    data <- data.frame(group,proportion,level)
    pos<-"fill"#"dodge"/"stack"/"fill"
    tmp_plot<-ggplot(data, aes(fill=group, y=level, x=proportion)) + ggtitle("Group proportions in each covariate level") +
      geom_bar(position=pos, stat="identity") + theme_bw() + scale_fill_viridis(discrete = T, option = "E", alpha=0.75)
    return(tmp_plot)
  }
}