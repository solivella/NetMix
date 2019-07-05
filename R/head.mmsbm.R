#' Identify nodes with most frequent membership in latent groups
#'
#' The function lists the nodes (optionally, node-time periods) that most frequently instantiate membership in each latent group.
#'
#' @param x An object of class \code{mmsbm}, a result of a call to \code{mmsbm}. 
#' @param n Numeric or integer; specifies how many units will be identified for each group.
#' @param t Optional vector of time periods to be used for assessing latent group membership.
#' @param node Logical; indicates whether latent group memberships should be averaged at the node level.  If FALSE, the function returns the node-time period units with highest estimated membership in each latent group.
#' @param t.correct Logical; indicates whether latent group memberships should be corrected for temporal trends.  If TRUE, the function returns the node-time period units with highest estimated membership in each latent group.  
#' @param ... Currently ignored
#'     
#' @return List of length \code{n.groups}.  Each entry contains a sorted vector of average latent membership probabilities of length \code{n}.
#
#' @method head mmsbm
#'
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 3 groups
#' set.seed(123)
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 3)
#' 
#' ## Show top 6 lawyers in each estimated latent block
#' head(lazega_mmsbm)


head.mmsbm <- function(x,  n=6, t=NULL, node=TRUE, t.correct=FALSE, ...){
  if(is.null(t)){t <- unique(x$monadic.data$`(tid)`)}
  Mem <- x$MixedMembership[,x$monadic.data[,"(tid)"] %in% t]
  if(!node){
    return(lapply(1:nrow(Mem), function(i){sort(Mem[i,], decreasing=TRUE)[1:n]}))
  }
  if(t.correct){
    Nodes <- x$monadic.data[,"(nid)"][x$monadic.data[,"(tid)"] %in% t]
    ts <- unlist(lapply(strsplit(colnames(Mem), "@"), "[[", 2))
    tm <- apply(Mem, 1, function(i){tapply(i, ts, mean)})
    Mem <- do.call(cbind, sapply(t, function(i){Mem[,ts==i] - tm[as.character(i),]}))
  }
  if(node){
    Nodes <- x$monadic.data[,"(nid)"][x$monadic.data[,"(tid)"] %in% t]
    node.mems <- t(do.call(cbind, lapply(unique(Nodes), function(i){
      rowMeans(as.matrix(Mem[,Nodes==i]))})))
    rownames(node.mems) <- as.character(unique(Nodes))
    return(lapply(1:x$n_blocks, function(i){
      node.mems[order(node.mems[,i], decreasing=T)[1:n],i]}))
  }
}

