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
#' @param node with \code{type == "membership"}; a character string specifying the node for which group membership should be plotted.
#' @param ... Currently ignored
#'
#' @return The requested plot object. 
#' 
#' @method plot mmsbm
#'
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aaylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 2 groups
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(seed = 123,
#'                                            hessian = FALSE))
#' 
#' ## Plot blockmodel as network
#' plot(lazega_mmsbm)
#' 



plot.mmsbm <- function(x, type="groups", FX=NULL, node=NULL, ...){ # network graph showing B-matrix
  if(type %in% c("blockmodel", "membership", "hmm")){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed to produce requested plot. Please install it.",
           call. = FALSE)
    }
  }
  
  if(type=="groups"){
    colRamp <- colorRamp(c("#DCDCDC","#808080","#000000"))
    g.mode <- ifelse(x$forms$directed, "directed", "undirected")
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
    v.size <- rowMeans(x$MixedMembership)*100 + 20
    radian.rescale <- function(x, start=0, direction=1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
      c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    loop.rads <- radian.rescale(x=1:x$n_blocks, direction=-1, start=0)
    loop.rads <- rep(loop.rads, times = times.arg)
    igraph::plot.igraph(block.G, main = "",
                        edge.width=4, edge.color=e.cols,  edge.curved = x$forms$directed, edge.arrow.size = 0.65,
                        edge.loop.angle = loop.rads,
                        vertex.size=v.size, vertex.color="white", vertex.frame.color="black",
                        vertex.label.font=2, vertex.label.cex=1, vertex.label.color="black",
                        layout = igraph::layout_in_circle)
    .bar.legend(colRamp, range(igraph::E(block.G)$weight))
  }
  
  if(type=="blockmodel"){
    x$dyadic.data$Y <- x$Y
    nodes <- unique(x$monadic.data$`(nid)`)
    MMat <- sapply(nodes,function(y){
      Dsub <- x$dyadic.data[x$dyadic.data$`(sid)`==y | x$dyadic.data$`(rid)`==y,]
      return(sapply(nodes, function(y){sum(Dsub$Y[Dsub$`(sid)`==y | Dsub$`(rid)`==y])}))
    })
    diag(MMat) <- 0
    clusters <- head(x, n=length(nodes))
    csort <- sort(sapply(nodes, function(y){which.max(sapply(clusters, "[[", y))}))
    corder <- unlist(sapply(1:x$n_blocks, function(z){sort(clusters[[z]][names(csort)[csort==z]], decreasing=T)}))
    MMat <- MMat[names(corder), names(corder)]
    plot(1, 1,
         xlim = c(.5, length(nodes) + .5),
         ylim = c(.5, length(nodes) + .5),
         main = "",
         xlab = "",
         ylab = "",
         type = "n", axes = FALSE)
    polygon.color <- c("white", "black")
    for (i in 1:length(nodes)) {
      for (t in 1:length(nodes)) {
        temp <- ifelse(MMat[i,t] > 0, 2, 1)
        polygon(c(.5 + t - 1, .5 + t, .5 + t, .5 + t - 1),
                length(nodes) - c(i-.5, i-.5, i+.5, i+.5),
                density = NA,
                border = polygon.color[temp],
                col = polygon.color[temp]) 
      } 
    }
    par(xpd=FALSE)
    for(i in 1:x$n_blocks){
      if(i < x$n_blocks){
        abline(h=length(nodes)-length(which(csort %in% 1:i))-.5, col="red", lty=2, lwd=2)
        abline(v=length(which(csort %in% 1:i))+.5, col="red", lty=2, lwd=2)
      }
    }
    
    v.size <- rowMeans(x$MixedMembership)
    bm <- x$BlockModel
    bm[upper.tri(bm)] <- NA
    dm <- data.frame(Sender = rep(paste("Group", 1:nrow(bm)), times = x$n_blocks),
                     Receiver = rep(paste("Group", 1:nrow(bm)), each = x$n_blocks),
                     Val = plogis(c(bm)))
    dm <- dm[complete.cases(dm),]
    dm$Sender  <- factor(dm$Sender, levels=rev(paste("Group", 1:nrow(bm))))
    p <- ggplot2::ggplot(ggplot2::aes_string(y = "Sender", x = "Receiver", fill="Val"), data = dm) +
      ggplot2::ggtitle("Edge Formation Between Blocs") + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::geom_tile(color = "white") + 
      ggplot2::theme_bw() +
      ggplot2::scale_size(guide='none') +
      ggplot2::scale_fill_gradient2(low = "#FEE0D2", mid = "#FB6A4A", high = "#99000D",
                                    midpoint = max(dm$Val)/2, limit = c(0,max(dm$Val)), name="Edge\nProbability")
    print(p)
  }
  
  if(type=="membership"){
    ifelse(is.null(node), 
           nr <- 1:nrow(x$monadic.data),
           nr <- which(x$monadic.data$`(nid)`==node))
    avgmems <- lapply(1:nrow(x$MixedMembership), function(y){
      tapply(x$MixedMembership[y,nr], x$monadic.data[nr,"(tid)"], mean)})
    avgmems <- as.data.frame(cbind(rep(unique(as.character(x$monadic.data[nr,"(tid)"])), nrow(x$MixedMembership)),unlist(avgmems),
                                   rep(1:nrow(x$MixedMembership), each=length(unique(x$monadic.data[nr,"(tid)"])))))
    colnames(avgmems) <- c("Time", "Avg.Membership", "Group")
    avgmems$Group <- factor(avgmems$Group, levels=length(unique(avgmems$Group)):1)
    if(class(avgmems$Avg.Membership) != "numeric"){avgmems$Avg.Membership <- as.numeric(as.character(avgmems$Avg.Membership))}
    if(class(avgmems$Time) != "numeric"){avgmems$Time <- as.numeric(as.character(avgmems$Time))}
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

