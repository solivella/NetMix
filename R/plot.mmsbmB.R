#' Plot output of mmsbm bipartite
#'
#' The function provides a variety of plotting options for a fitted \code{mmsbmB} object.
#'
#' @param fm An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}.
#' @param type character string denoting the type of plot. The default, "\code{blockmodel}," plots the estimated matrix of 
#' group by group edge formation probabilities as a network graph.  "\code{membership}" plots average membership in
#' each latent group by time period. "\code{effect}" provides a series of plots showing the estimated effect 
#' of a shfit in monadic covariate values.
#' @param FX with type = "effect"; a list resulting from a call to \code{covFX}.
#'


plot.mmsbmB <- function(x, type="groups", FX=NULL,...){ # network graph showing B-matrix
  if(type %in% c("blockmodel", "membership", "hmm")){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed to produce requested plot. Please install it.",
           call. = FALSE)
    }
  }
  
  if(type=="groups"){
    colRamp <- colorRamp(c("#DCDCDC","#808080","#000000"))
    g.mode<-"undirected"
    adj_mat <- x$BlockModel
    dimnames(adj_mat) <- list(paste("G",1:nrow(adj_mat), sep=""),
                              paste("H", 1:ncol(adj_mat), sep=""))  
    block.G <- igraph::graph.incidence(plogis(adj_mat), weighted=TRUE)
    
    e.weight <- (1/diff(range(igraph::E(block.G)$weight))) * (igraph::E(block.G)$weight - max(igraph::E(block.G)$weight)) + 1
    e.cols <- rgb(colRamp(e.weight), maxColorValue = 255)
    v.size <- c(rowMeans(x$`MixedMembership 1`)*100 + 25, rowMeans(x$`MixedMembership 2`)*100 + 25)
    igraph::plot.igraph(block.G, main = "",
                        edge.width=4, edge.color=e.cols,  edge.curved = F, 
                        edge.arrow.size = 0.65,
                        vertex.size=v.size, 
                        vertex.color="white", vertex.frame.color="black",
                        vertex.label.font=2, vertex.label.cex=1, vertex.label.color="black",
                        layout = igraph::layout.bipartite)
    .bar.legend(colRamp, range(igraph::E(block.G)$weight))
  }
  
  if(type=="membership"){
    avgmems <- lapply(1:nrow(x$MixedMembership), function(x){
      tapply(x$MixedMembership[x,], x$monadic.data[,"(tid)"], mean)})
    avgmems <- as.data.frame(cbind(rep(unique(as.character(x$monadic.data[,"(tid)"])), nrow(x$MixedMembership)),unlist(avgmems),
                                   rep(1:nrow(x$MixedMembership), each=length(unique(x$monadic.data[,"(tid)"])))))
    colnames(avgmems) <- c("Time", "Avg.Membership", "Group")
    avgmems$Group <- factor(avgmems$Group, levels=length(unique(avgmems$Group)):1)
    if(class(avgmems$Avg.Membership) == "factor"){avgmems$Avg.Membership <- as.numeric(as.character(avgmems$Avg.Membership))}
    if(class(avgmems$Time) == "factor"){avgmems$Time <- as.numeric(as.character(avgmems$Time))}
    return(ggplot2::ggplot() + 
             ggplot2::geom_area(ggplot2::aes_string(y = "Avg.Membership", x = "Time", fill="Group"), data = avgmems,
                                stat="identity", position="stack") + 
             ggplot2::guides(fill=ggplot2::guide_legend(title="Group")))
  }
  
  if(type=="effect"){
    stopifnot(is.list(FX))
    cov <- strsplit(names(FX)[1], " ")[[1]][5]
    ymax <- max(hist(FX[[5]])[["counts"]])
    hist(FX[[5]], main=paste("Distribution of Marginal Effects:", strsplit(names(FX)[1], " ")[[1]][5]),
         xlab=paste("Effect of", cov, "on Pr(Edge Formation)"))
    
    plot(unique(x$dyadic.data[,"(tid)"]), tapply(FX[[5]], x$dyadic.data[,"(tid)"], mean), type="o",
         xlab="Time", ylab=paste("Effect of", cov, "on Pr(Edge Formation)"), main="Marginal Effect over Time")
    
    nodenames <- names(sort(table(x$monadic.data[,"(nid)"]), decreasing=TRUE))
    nodes <- sort(FX[[3]])[names(sort(FX[[3]])) %in% nodenames]
    plot(1, type="n", xlab="Node-Level Estimated Effect", ylab="", 
         xlim=c(min(nodes), max(nodes) + 0.001),
         ylim = c(0, length(nodes)), yaxt="n")
    for(i in 1:length(nodes)){
      points(nodes[i],i, pch=19)
      text(nodes[i],i, names(nodes)[i], pos=4, cex=0.7)
    }
  }
  
  if(type=="hmm"){
    hms <- as.data.frame(do.call(rbind, lapply(1:nrow(x$Kappa), function(x){
      cbind(1:ncol(x$Kappa), x$Kappa[x,], x)
    })))
    colnames(hms) <- c("Time", "Kappa", "State")
    hms$State <- as.factor(hms$State)
    return(ggplot2::ggplot() + 
             ggplot2::geom_area(ggplot2::aes_string(y = "Kappa", x = "Time", fill="State"), data = hms,
                                stat="identity", position="stack") + 
             ggplot2::guides(fill=ggplot2::guide_legend(title="HMM State")))
  }
  
  if(type=="block"){
    adj_mat <- x$BlockModel
    dimnames(adj_mat) <- list(paste("G",1:nrow(adj_mat), sep=""),
                              paste("H", 1:ncol(adj_mat), sep=""))  
    melt_block<-melt(plogis(adj_mat))
    colnames(melt_block)<-c("G","H","Probability")
    return(ggplot2::ggplot(data = melt_block, aes(x=H, y=G, fill=`Probability`)) + 
             ggplot2::geom_tile() + ggplot2::geom_text(aes(label=round(`Probability`,3)), col="white"))
  }
}