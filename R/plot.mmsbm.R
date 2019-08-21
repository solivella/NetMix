#' Various visualization tools for 'mmsbm' objects
#'
#' The function provides a variety of plotting options for a fitted \code{mmsbm} object.
#'
#' @param x An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param type character string denoting the type of plot. The default, "\code{groups}," plots the estimated matrix of 
#' group by group edge formation probabilities as a network plot, with nodes representing groups (sized proportional to relative membership) 
#' and edge colors encoding probability of between-group ties. "\code{blockmodel}" plots the same information,
#' but using a tile plot instead of a network plot.  "\code{membership}" plots average membership in
#' each latent group by time period. "\code{effect}" provides a series of plots showing the estimated effect 
#' of a shfit in monadic covariate values.
#' @param FX with \code{type == "effect"}; a list resulting from a call to \code{covFX}.
#' @param ... Currently ignored
#'
#' @return The requested plot object. 
#' 
#' @method plot mmsbm
#'
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (adelinel@@princeton.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 2 groups
#' set.seed(123)
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2)
#' 
#' ## Plot blockmodel as network
#' plot(lazega_mmsbm)
#' 



plot.mmsbm <- function(x, type="groups", FX=NULL, ...){ # network graph showing B-matrix
  if(type %in% c("blockmodel", "membership", "hmm")){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed to produce requested plot. Please install it.",
           call. = FALSE)
    }
  }
  
  if(type=="groups"){
    colRamp <- colorRamp(c("#DCDCDC","#808080","#000000"))
    g.mode <- ifelse(eval(x$forms$directed), "directed", "undirected")
    adj_mat <- x$BlockModel
    dimnames(adj_mat) <- list(paste("G",1:nrow(adj_mat), sep=""),
                              paste("G", 1:ncol(adj_mat), sep=""))
    block.G <- igraph::graph.adjacency(plogis(adj_mat), mode=g.mode, weighted=TRUE)
    e.weight <- (1/diff(range(igraph::E(block.G)$weight))) * (igraph::E(block.G)$weight - max(igraph::E(block.G)$weight)) + 1
    e.cols <- rgb(colRamp(e.weight), maxColorValue = 255)
    times.arg <- if(g.mode == "directed") {
      x$n_blocks
    } else {
      rev(seq_len(x$n_blocks))
    }
    v.size <- rowMeans(x$MixedMembership)*100 + 30
    loop.rads <- rep(seq(0, -1.75*pi, length.out = x$n_blocks),
                     times = times.arg)
    igraph::plot.igraph(block.G, main = "",
                        edge.width=4, edge.color=e.cols,  edge.curved = x$directed, edge.arrow.size = 0.65,
                        edge.loop.angle = loop.rads,
                        vertex.size=v.size, vertex.color="gray80", vertex.frame.color="black",
                        vertex.label.font=2, vertex.label.cex=1, vertex.label.color="black",
                        layout = igraph::layout_in_circle)
    .bar.legend(colRamp)
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
}

