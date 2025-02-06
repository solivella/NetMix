#' Plot output of mmsbm bipartite
#'
#' The function provides a variety of plotting options for a fitted \code{mmsbmB} object.
#'
#' @param fm An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}.
#' @param type character string denoting the type of plot. The default, "\code{blockmodel}," plots the estimated matrix of 
#' group by group edge formation probabilities as a network graph.  "\code{membership}" plots average membership in
#' each latent group by time period. "\code{effect}" provides a series of plots showing the estimated effect 
#' of a shfit in monadic covariate values.
#' @param FX with type = "effect"; a list resulting from a call to \code{covFXB}.
#' @param family with type = "effect"; integer 1 or 2 for whether the effect is meant for Family 1 or Family 2 nodes.
#' @param nodelabel with type = "effect"; list of node names for node-effect plot, 
#' must be same length as number of nodes and in the original order of nid passed to mmsbmB.
#' @param blocklabel with type = "blockmodel"; list of block labels, with one character vector per family.


plot.mmsbmB <- function(x, type="groups", FX=NULL, family=1, nodelabel=NULL,...){ # network graph showing B-matrix
  if(type %in% c("groups", "membership", "hmm", "block")){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed to produce requested plot. Please install it.",
           call. = FALSE)
    }
  }
  if(type %in% c("groups")){
    if (!requireNamespace("ggnetwork", quietly = TRUE)) {
      stop("Package \"ggnetwork\" needed to produce requested plot. Please install it.",
           call. = FALSE)
    }
  }
  
  all_args <- list(...)
  if(type=="groups"){
    require(igraph, quietly = TRUE)
    require(ggraph, quietly = TRUE)
    adj_mat <- x$BlockModel
    n_vertex <- ifelse(x$bipartite, sum(dim(adj_mat)), nrow(adj_mat))
    n_edge <- prod(dim(adj_mat))
    if(is.null(all_args$vertex.label)){
      col_prefix <- ifelse(x$bipartite, "H", "G")
      dimnames(adj_mat) <- list(paste("G",1:nrow(adj_mat), sep=""),
                                paste(col_prefix,1:nrow(adj_mat), sep=""))
      vertex.label <- rownames(adj_mat)
    } else {
      vertex.label <- all_args$vertex.label
    }
    if(is.null(all_args$vertex.color)){
      vertex.color <- rep("gray50", n_vertex)
    } else {
      vertex.color <- all_args$vertex.color
    }
    if(x$bipartite){
      sizes <- c(rowMeans(x$`MixedMembership1`)*100, rowMeans(x$`MixedMembership2`)*100)
      dir <- 45
      layout_name <- "bipartite"
    } else {
      sizes <- rowMeans(x$`MixedMembership1`)*100
      dir <- seq(45, 360, length.out = n_edge)
      layout_name <- "circle"
    }
    graph_fun <- ifelse(x$bipartite, igraph::graph_from_biadjacency_matrix, igraph::graph_from_adjacency_matrix)
    block.G <- graph_fun(plogis(adj_mat), weighted=TRUE) %>% 
      set_vertex_attr("v.lab",value = vertex.label) %>% 
      set_vertex_attr("v.col",value = vertex.color) %>% 
      set_vertex_attr("MM",value = sizes)
    
    
    bm_lo <- create_layout(block.G,  layout = "igraph", algorithm=layout_name)
    if(x$bipartite){
      adj_x <- 0
      adj_y <- ifelse(V(block.G)$type, -0.17, 0.17)
    } else {
      adj_x <- bm_lo$x *0.25
      adj_y <- bm_lo$y *0.25
    }
    return(ggraph(bm_lo) +
             geom_edge_link(aes(color = weight), linewidth=1.5) +
             geom_edge_loop(aes(color = weight,
                                direction = dir,
                                span = 60,
                                strength=0.4),
                            linewidth=1.5) +
             geom_node_point(aes(size=MM, fill=v.col, color=v.col),
                             show.legend = FALSE) +
             scale_edge_color_gradient("Edge\nProbability",
                                       low = "gray90", high = "black", 
                                       limits=c(0,1)) +
             scale_size_area(max_size = 15, guide="none") +
             geom_node_text(aes(label = v.lab),
                            fontface = "bold",
                            size = 5,
                            nudge_x = adj_x,
                            nudge_y = adj_y) +
             scale_fill_manual(values = vertex.color) + 
             scale_colour_manual(values = vertex.color) + 
             theme_void() +
             theme(legend.justification = ifelse(x$bipartite, "center","bottom"),
                   legend.title = element_text(size=12)) +
             coord_cartesian(clip="off"))
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
    #hist(FX[[5]], main=paste("Distribution of Marginal Effects:", strsplit(names(FX)[1], " ")[[1]][5]),
    #xlab=paste("Effect of", cov, "on Pr(Edge Formation)"))
    #plot(unique(x$dyadic.data[,x$forms$timeID]), tapply(FX[[5]], x$dyadic.data[,x$forms$timeID], mean), type="o",
    #xlab="Time", ylab=paste("Effect of", cov, "on Pr(Edge Formation)"), main="Marginal Effect over Time")
    
    if(family==1){ 
      monadic.data<- x$monadic.data[[1]]
      nid <- x$forms$senderID
      if(!nid%in%names(monadic.data)){nid<-"(nid1)"}
    } else{ 
      monadic.data <- x$monadic.data[[2]]
      nid <- x$forms$receiverID
      if(!nid%in%names(monadic.data)){nid<-"(nid2)"}
    }
    
    ##set up alternative labels for nodes -- these need to be in the same starting order as nid in monadic.data
    if(!is.null(nodelabel)){
      tmp_labels<-data.frame(nodelabel=nodelabel,nodenames=monadic.data[,nid])
    }
    nodenames <- names(sort(table(monadic.data[,nid]), decreasing=TRUE))
    nodes <- sort(FX[[3]])[names(sort(FX[[3]])) %in% nodenames]
    if(!is.null(nodelabel)){
      nodes_labels<- tmp_labels$nodelabel[match(names(nodes),tmp_labels$nodenames)]
    }else{
      nodes_labels<-names(nodes)
    }
    plot(1, type="n", xlab="Node-Level Estimated Effect", ylab="", 
         xlim=c(min(nodes), max(nodes) + 0.001),
         ylim = c(0, length(nodes)), yaxt="n")
    for(i in 1:length(nodes)){
      points(nodes[i],i, pch=19)
      text(nodes[i],i, nodes_labels[i], pos=4, cex=0.5)
    }
  }
  
  if(type=="effectchange"){
    stopifnot(is.list(FX))
    cov <- strsplit(names(FX)[1], " ")[[1]][5]
    #ymax <- max(hist(FX[[5]])[["counts"]])
    #hist(FX[[5]], main=paste("Distribution of Marginal Effects:", strsplit(names(FX)[1], " ")[[1]][5]),
    #xlab=paste("Effect of", cov, "on Pr(Edge Formation)"))
    #plot(unique(x$dyadic.data[,x$forms$timeID]), tapply(FX[[5]], x$dyadic.data[,x$forms$timeID], mean), type="o",
    #xlab="Time", ylab=paste("Effect of", cov, "on Pr(Edge Formation)"), main="Marginal Effect over Time")
    
    if(family==1){ 
      monadic.data<- x$monadic.data[[1]]
      nid <- x$forms$senderID
      if(!nid%in%names(monadic.data)){nid<-"(nid1)"}
    } else{ 
      monadic.data <- x$monadic.data[[2]]
      nid <- x$forms$receiverID
      if(!nid%in%names(monadic.data)){nid<-"(nid2)"}
    }
    ##set up alternative labels for nodes -- these need to be in the same starting order as nid in monadic.data
    if(!is.null(nodelabel)){
      tmp_labels<-data.frame(nodelabel=nodelabel,nodenames=monadic.data[,nid])
    }
    nodenames <- names(sort(table(monadic.data[,nid]), decreasing=TRUE))
    tmp_order<- names(sort(FX[[3]]))
    predicted_nodes <- FX[[6]][match(tmp_order,names(FX[[6]]))]
    orig_nodes <- FX[[7]][match(tmp_order,names(FX[[7]]))]
    if(!is.null(nodelabel)){
      nodes_labels<- tmp_labels$nodelabel[match(tmp_order,tmp_labels$nodenames)]
    }else{
      nodes_labels<-tmp_order
    }
    plot(1, type="n", xlab="Node-Level Estimated Change", ylab="", 
         xlim=c(min(c(predicted_nodes,orig_nodes)), max(c(predicted_nodes,orig_nodes)) + 0.001),
         ylim = c(0, length(predicted_nodes)), yaxt="n")
    for(i in 1:length(predicted_nodes)){
      points(predicted_nodes[i],i, pch=19)
      points(orig_nodes[i],i, pch=19, col="gray")
      lines(x=c(predicted_nodes[i],orig_nodes[i]),y=rep(i,2),lty='dashed',col="gray")
      text(max(c(predicted_nodes[i],orig_nodes[i])),i, nodes_labels[i], pos=4, cex=0.5)
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
             ggplot2::geom_tile(aes(fill = `Probability`)) + ggplot2::geom_text(aes(label=round(`Probability`,3)), col="dodgerblue4") +
             ggplot2::scale_fill_gradient(low = "white", high = "black"))
  }
}
