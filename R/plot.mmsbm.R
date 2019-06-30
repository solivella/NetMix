#' Various visualization tools for 'mmsbm' objects
#'
#' The function provides a variety of plotting options for a fitted \code{mmsbm} object.
#'
#' @param x An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param type character string denoting the type of plot. The default, "\code{blockmodel}," plots the estimated matrix of 
#' group by group edge formation probabilities as a tile plot, with tiles proportional to group sizes.  "\code{membership}" plots average membership in
#' each latent group by time period. "\code{effect}" provides a series of plots showing the estimated effect 
#' of a shfit in monadic covariate values.
#' @param FX with \code{type == "effect"}; a list resulting from a call to \code{covFX}.
#' @param ... Currently ignored
#'
#' @return A plot object. 
#' 
#' @method plot mmsbm
#'
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' 


plot.mmsbm <- function(x, type="blockmodel", FX=NULL, ...){ # network graph showing B-matrix
  if(type=="blockmodel"){
    mode <- ifelse(eval(x$forms$directed), "directed", "undirected")
    v.size <- rowMeans(x$MixedMembership)
    
    bm <- x$BlockModel
    if(mode!="directed") {
      bm[upper.tri(bm)] <- NA
    }
    
    dm <- data.frame(Sender = rep(rownames(x$BlockModel), times = x$n_blocks),
                     Receiver = rep(colnames(x$BlockModel), each = x$n_blocks),
                     Val = plogis(c(bm)),
                     Height = rep(v.size*(x$n_blocks-0.35), x$n_blocks),
                     Width = rep(v.size*(x$n_blocks-0.35), each=x$n_blocks))
    
    dm <- dm[complete.cases(dm),]
    
    dm$Sender  <- factor(dm$Sender, levels=rev(rownames(x$BlockModel)))
    
    ggplot2::ggplot(aes_string(y = "Sender", x ="Receiver", fill="Val", width="Width", height="Height"),
                    data = dm) + 
      ggplot2::geom_tile(color = "white") + theme_bw()+
      ggplot2::scale_size(guide='none') +
      ggplot2::scale_fill_gradient2(low = "#deebf7", mid = "#9ecae1", high = "#3182bd",
                           midpoint = 0.5, limit = c(0,1), name="Edge\nProbability")
  }
  
  if(type=="groups"){
    colRamp <- colorRamp(c("#ffffcc","#fd8d3c","#800026"))
    mode <- ifelse(eval(x$forms$directed), "directed", "undirected")
    adj_mat <- x$BlockModel
    dimnames(adj_mat) <- list(paste("G",1:nrow(adj_mat), sep=""),
                              paste("G", 1:ncol(adj_mat), sep=""))
    block.G <- igraph::graph.adjacency(plogis(adj_mat), mode=mode, weighted=TRUE)
    e.weight <- igraph::E(block.G)$weight
    e.cols <- rgb(colRamp(e.weight), maxColorValue = 255)
    v.size <- rowMeans(x$MixedMembership)*100
    layout(matrix(1:2,ncol=2), widths = c(2,1), heights = c(1,1))
    opar <- par(mar=c(0,0,0,0)+.75)
    igraph::plot.igraph(block.G, main = "",
         edge.width=4, edge.color=e.cols,  edge.curved = 0.5, edge.arrow.size = .65,
         vertex.size=v.size, vertex.color="gray80", vertex.frame.color="black",
         vertex.label.font=2, vertex.label.cex=1,
         layout = igraph::layout_nicely)
    par(opar)
    legend_image <- as.raster(matrix(rgb(colRamp(seq(1,0,length.out=50)), maxColorValue = 255), ncol=1))
    plot(c(0,2.5),c(0.3,.7),type = 'n', axes = FALSE ,xlab = '', ylab = '')
    title('Edge\nprobability', cex.main=0.9, adj=0, line=-0.005)
    text(x=1.9, y = seq(0.3,.7,length.out = 5), labels = seq(0,1,length.out =5), cex = 0.75)
    rasterImage(legend_image, 0, 0.3, 1,0.7)
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
    ggplot2::ggplot() + 
      ggplot2::geom_area(aes_string(y = "Avg.Membership", x = "Time", fill="Group"), data = avgmems,
                stat="identity", position="stack") + 
      ggplot2::guides(fill=guide_legend(title="Group"))
  }
  
  if(type=="effect"){
    stopifnot(is.list(FX))
    cov <- strsplit(names(FX)[1], " ")[[1]][5]
    ymax <- max(hist(FX[[5]])[["counts"]])
    hist(FX[[5]], main=paste("Distribution of Marginal Effects:", strsplit(names(FX)[1], " ")[[1]][5]),
         xlab=paste("Effect of", cov, "on Pr(Edge Formation)"))
    
    plot(unique(x$dyadic.data[,"(tid)"]), tapply(FX[[5]], x$dyadic.data[,"(tid)"], mean), type="o",
         xlab="Time", ylab=paste("Effect of", cov, "on Pr(Edge Formation)"), main="Marginal Effect over Time")
    
    nodenames <- names(sort(table(x$monadic.data[,"(nid)"]), decreasing=T))
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
    ggplot2::ggplot() + 
      ggplot2::geom_area(aes_string(y = "Kappa", x = "Time", fill="State"), data = hms,
                stat="identity", position="stack") + 
      ggplot2::guides(fill=guide_legend(title="HMM State"))
  }
}

