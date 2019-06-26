#' Plot output of dynMMSBM
#'
#' The function provides a variety of plotting options for a fitted \code{mmsbm} object.
#'
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param type character string denoting the type of plot. The default, "\code{blockmodel}," plots the estimated matrix of 
#' group by group edge formation probabilities as a tile plot, with tiles proportional to group sizes.  "\code{membership}" plots average membership in
#' each latent group by time period. "\code{effect}" provides a series of plots showing the estimated effect 
#' of a shfit in monadic covariate values.
#' @param FX with \code{type == "effect"}; a list resulting from a call to \code{covFX}.
#'


plot.mmsbm <- function(fm, type="blockmodel", FX=NULL){ # network graph showing B-matrix
  if(type=="blockmodel"){
    mode <- ifelse(eval(fm$forms$directed), "directed", "undirected")
    # require("igraph", quietly=TRUE)
    # block.G <- graph.adjacency(exp(fm$BlockModel) / (1 + exp(fm$BlockModel)), mode=mode, weighted=TRUE)
    # e.weight <- E(block.G)$weight*100
    # if(any(e.weight < 0.1)){
    #   e.weight <- e.weight*2
    # }
    # if(any(e.weight > 18)){
    #   e.weight[e.weight > 18] <- 18}
    v.size <- rowMeans(fm$MixedMembership)
    # plot(block.G, main = "Edge Formation across Clusters",
    #     edge.width=e.weight, vertex.size=v.size,
    #     layout = layout_in_circle)
    
    bm <- fm$BlockModel
    if(mode!="directed") {
      bm[upper.tri(bm)] <- NA
    }
    
    dm <- data.frame(Sender = rep(rownames(fm$BlockModel), times = fm$n_blocks),
                     Receiver = rep(colnames(fm$BlockModel), each = fm$n_blocks),
                     Val = plogis(c(bm)),
                     Height = rep(v.size*(fm$n_blocks-0.35), fm$n_blocks),
                     Width = rep(v.size*(fm$n_blocks-0.35), each=fm$n_blocks))
    
    dm <- dm[complete.cases(dm),]
    
    dm$Sender  <- factor(dm$Sender, levels=rev(rownames(fm$BlockModel)))
    
    require(ggplot2)
    p <- ggplot(aes(y = Sender, x = Receiver, fill=Val, width=Width, height=Height), data = dm) + ggtitle("Edge Formation Between Blocs") + theme(plot.title = element_text(hjust = 0.5)) +
      geom_tile(color = "white") + theme_bw()+
      scale_size(guide='none') +
      scale_fill_gradient2(low = "#deebf7", mid = "#9ecae1", high = "#3182bd",
                           midpoint = 0.5, limit = c(0,1), name="Edge\nProbability")
    print(p)
  }
  
  if(type=="groups"){
    require("igraph", quietly=TRUE)
    mode <- ifelse(eval(fm$forms$directed), "directed", "undirected")
    block.G <- graph.adjacency(exp(fm$BlockModel) / (1 + exp(fm$BlockModel)), mode=mode, weighted=TRUE)
    e.weight <- E(block.G)$weight*100
    if(any(e.weight < 0.1)){
      e.weight <- e.weight*2
    }
    if(any(e.weight > 18)){
      e.weight[e.weight > 18] <- 18}
    v.size <- rowMeans(fm$MixedMembership)*100
    plot(block.G, main = "Edge Formation Between Latent Groups",
         edge.width=e.weight, vertex.size=v.size,
         vertex.label.font=2, vertex.label.cex=1.1,
         layout = layout_with_gem)
  }
  
  if(type=="membership"){
    avgmems <- lapply(1:nrow(fm$MixedMembership), function(x){
      tapply(fm$MixedMembership[x,], fm$monadic.data[,"(tid)"], mean)})
    avgmems <- as.data.frame(cbind(rep(unique(as.character(fm$monadic.data[,"(tid)"])), nrow(fm$MixedMembership)),unlist(avgmems),
                                   rep(1:nrow(fm$MixedMembership), each=length(unique(fm$monadic.data[,"(tid)"])))))
    colnames(avgmems) <- c("Time", "Avg.Membership", "Group")
    avgmems$Group <- factor(avgmems$Group, levels=length(unique(avgmems$Group)):1)
    if(class(avgmems$Avg.Membership) == "factor"){avgmems$Avg.Membership <- as.numeric(as.character(avgmems$Avg.Membership))}
    if(class(avgmems$Time) == "factor"){avgmems$Time <- as.numeric(as.character(avgmems$Time))}
    require(ggplot2)
    p <- ggplot() + ggtitle("Average Group Membership Over Time") + theme(plot.title = element_text(hjust = 0.5)) +
          geom_area(aes(y = Avg.Membership, x = Time, fill=Group), data = avgmems,
                stat="identity", position="stack")  + guides(fill=guide_legend(title="Group"))
    print(p)
  }
  
  if(type=="effect"){
    stopifnot(is.list(FX))
    cov <- strsplit(names(FX)[1], " ")[[1]][5]
    ymax <- max(hist(FX[[5]])[["counts"]])
    hist(FX[[5]], main=paste("Distribution of Marginal Effects:", strsplit(names(FX)[1], " ")[[1]][5]),
         xlab=paste("Effect of", cov, "on Pr(Edge Formation)"))
    
    plot(unique(fm$dyadic.data[,"(tid)"]), tapply(FX[[5]], fm$dyadic.data[,"(tid)"], mean), type="o",
         xlab="Time", ylab=paste("Effect of", cov, "on Pr(Edge Formation)"), main="Marginal Effect over Time")
    
    nodenames <- names(sort(table(fm$monadic.data[,"(nid)"]), decreasing=T))
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
    hms <- as.data.frame(do.call(rbind, lapply(1:nrow(fm$Kappa), function(x){
      cbind(1:ncol(fm$Kappa), fm$Kappa[x,], x)
      })))
    colnames(hms) <- c("Time", "Kappa", "State")
    require(ggplot2)
    p <- ggplot() + ggtitle("Hidden Markov State") + theme(plot.title = element_text(hjust = 0.5)) +
      geom_area(aes(y = Kappa, x = Time, fill=factor(State)), data = hms,
                stat="identity", position="stack")  + guides(fill=guide_legend(title="State"))
    print(p)
  }
}

