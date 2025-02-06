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



plot.mmsbm <- function(x, type="groups", FX=NULL, ...){ # network graph showing B-matrix
  if(type %in% c("blockmodel", "membership", "hmm")){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package \"ggplot2\" needed to produce requested plot. Please install it.",
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
      set_vertex_attr("MM",value = sizes) %>% 
      set_edge_attr("MM",value = rep(sizes, each=n_vertex))
    
    
    bm_lo <- create_layout(block.G,  layout = "igraph", algorithm=layout_name)
    if(x$bipartite){
      adj_x <- 0
      adj_y <- ifelse(V(block.G)$type, -0.1, 0.1)
    } else {
      adj_x <- bm_lo$x *0.25
      adj_y <- bm_lo$y *0.25
    }
    return(ggraph(bm_lo) +
             geom_edge_link(aes(color = weight), linewidth=1.5) +
             geom_edge_loop(aes(color = weight,
                                direction = dir,
                                span = 60,
                                strength=1),
                            linewidth=1.5) +
             geom_node_point(aes(size=MM, fill=v.lab, color=v.lab),
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
             theme(legend.justification = ifelse(x$bipartite, "center","top"),
                   legend.title = element_text(size=12)) +
             coord_cartesian(clip="off"))
  }
  
  if(type=="membership"){
    avgmems <- lapply(1:nrow(x$MixedMembership1), function(x){
      tapply(x$MixedMembership1[x,], x$monadic.data[[1]][,"(tid)"], mean)})
    avgmems <- as.data.frame(cbind(rep(unique(as.character(x$monadic.data[[1]][,"(tid)"])), nrow(x$MixedMembership1)),unlist(avgmems),
                                   rep(1:nrow(x$MixedMembership1), each=length(unique(x$monadic.data[[1]][,"(tid)"])))))
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
    
    nodenames <- names(sort(table(x$monadic.data[[1]][,"(nid)"]), decreasing=TRUE))
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
