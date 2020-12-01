#' Output top members of each latent group based on estimated mmsbm model
#'
#' The function produces a data.frame with rows=group*n, cols=group, names, probas
#'
#' @param membership Item `MixedMembership 1` or `MixedMembership 2` from object of class \code{mmsbmB}.
#' @param n Number of top members to return \code{integer} object. 
#' @param names Vector of names for nodes \code{vector} object. 
#'     
#' @author Kosuke Imai (imai@@harvard.edu), Adeline Lo (aylo@@wisc.edu), Santiago Olivella (olivella@@unc.edu), Tyler Pratt (tyler.pratt@@yale.edu)
#' 
#' @example tests/Examples/cosponsorhip.rmd
#' 
#' 
topmmB<-function(membership,n,names){
  group<-nrow(membership)
  res<-vector("list",length=group)
  for(i in 1:group){
    g<-membership[i,]
    names(g)<-names
    g<-sort(g,decreasing = TRUE)[1:n]
    res[[i]]<-data.frame(group=i,names=names(g),probability=g)
  }
  res<-do.call(rbind,res) #returns data.frame with rows=group*n, cols=group, names, probas
  return(res)
}
