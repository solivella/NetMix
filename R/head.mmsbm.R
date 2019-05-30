#' Identify nodes with most frequent membership in latent groups
#'
#' The function lists the nodes (optionally, node-time periods) that most frequently instantiate membership in each latent group.
#'
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}. 
#' @param n numeric or integer; specifies how many units will be identified for each group.
#' @param t optional vector of time periods to be used for assessing latent group membership.
#' @param node logical; indicates whether latent group memberships should be averaged at the node level.  If FALSE, the function returns the node-time period units with highest estimated membership in each latent group.
#' @param t.correct logical; indicates whether latent group memberships should be corrected for temporal trends.  If TRUE, the function returns the node-time period units with highest estimated membership in each latent group.  
#'
#'     
#' @return List of length \code{n.groups}.  Each entry contains a sorted vector of average latent membership probabilities of length \code{n}.
#
#' 



head.mmsbm <- function(fm,  n=10, t=NULL, node=TRUE, t.correct=FALSE){
  if(is.null(t)){t <- unique(fm$monadic.data$`(tid)`)}
  Mem <- fm$MixedMembership[,fm$monadic.data[,"(tid)"] %in% t]
  if(!node){
    return(lapply(1:nrow(Mem), function(x){sort(Mem[x,], decreasing=T)[1:n]}))
  }
  if(t.correct){
    Nodes <- fm$monadic.data[,"(nid)"][fm$monadic.data[,"(tid)"] %in% t]
    ts <- unlist(lapply(strsplit(colnames(Mem), "@"), "[[", 2))
    tm <- apply(Mem, 1, function(x){tapply(x, ts, mean)})
    Mem <- do.call(cbind, sapply(t, function(x){Mem[,ts==x] - tm[as.character(x),]}))
  }
  if(node){
    Nodes <- fm$monadic.data[,"(nid)"][fm$monadic.data[,"(tid)"] %in% t]
    node.mems <- t(do.call(cbind, lapply(unique(Nodes), function(x){
      rowMeans(as.matrix(Mem[,Nodes==x]))})))
    rownames(node.mems) <- as.character(unique(Nodes))
    return(lapply(1:fm$n_blocks, function(x){
      node.mems[order(node.mems[,x], decreasing=T)[1:n],x]}))
  }
}

