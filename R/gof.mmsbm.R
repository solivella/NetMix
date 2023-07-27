#' Posterior predictive checks using structural network characteristics
#'
#' The function generates a variety of plots that serve as posterior predictive checks on the goodness of fit of a fitted \code{mmsbm} object.
#'
#' Goodness of fit of network models has typically been established by evaluating how the structural characteristics of predicted networks 
#' compare to those of the observed network. When estimated in a Bayesian framework, this approach is equivalent to 
#' conducting posterior predictive checks on these structural quantities of interest. When \code{new.data.dyad} and/or \code{new.data.monad} are
#' passed that are different from those used in estimation, this is equivalent to conducting posterior predictive checks out-of-sample. 
#' 
#' The set of structural features used to determine goodness of fit is somewhat arbitrary, and chosen mostly to incorporate various 
#' first order, second order, and (to the extent possible) third-order characteristics of the network. "Geodesics" focuses on the distribution over 
#' observed and predicted geodesic distances between nodes; "Indegree" and "Outdegree" focuses on the distribution over incoming and outgoing connections 
#' per node; "3-motifs" focus on a distribution over possible connectivity patterns between triads (i.e., the triadic census); "Dyad Shared Partners" focuses on the distribution
#' over the number of shared partners between any two dyads; "Edge Shared Partners" is similarly defined, but w.r.t. edges, rather than dyads; and finally
#' "Incoming K-stars" focuses on a frequency distribution over stars with k=1,... spokes. 
#' 
#' Obtaining samples of the last three structural features can be very computationally expensive, and is discouraged on networks with more than 50
#' nodes.
#' 
#' For bipartite networks, first order structural features can be defined by family. Accordingly, \code{gof_stat} also accepts "Bipartite Degree" as
#' an input when \code{x} is of class \code{mmsbmBipartite}, which will generate a separate gof plot for each family.   
#'  
#'  
#' @param x An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param gof_stat Character vector. Accepts any subset from "Geodesics","Degree", "Indegree", "Outdegree", "3-Motifs",
#'                 "Dyad Shared Partners", "Edge Shared Partners", and "Incoming K-stars". See details.
#' @param level Double. Level of credible interval for posterior predictive distribution around structural quantities of interest. 
#' @param samples Integer. Number of sampled networks from model's posterior predictive using \code{\link{simulate.mmsbm}}.
#' @param new.data.dyad See \code{\link{simulate.mmsbm}}. Enables out-of-sample checking.
#' @param new.data.monad See \code{\link{simulate.mmsbm}}. Enables out-of-sample checking.
#' @param seed See \code{\link{simulate.mmsbm}}.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot} object.
#'    
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' 
#' ## Estimate model with 2 groups
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(seed = 123,
#'                                            hessian = FALSE))
#' 
#' ## Plot observed (red) and simulated (gray) distributions over 
#' ## indegrees
#' ## (typically a larger number of samples would be taken) 
#' gof(lazega_mmsbm, gof_stat = "Indegree", samples = 2)
#'

gof <- function (x, ...) {
  UseMethod("gof", x)
}

#' @method gof mmsbm
#' @rdname gof
gof.mmsbm <- function(x,
                      gof_stat = c("Geodesics","Degree"),
                      level = 0.95,
                      samples = 50,
                      new.data.dyad = NULL,
                      new.data.monad  = NULL, 
                      seed = NULL,
                      ...
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed to produce plot of GOF statistics. Please install it.",
         call. = FALSE)
  }
  
  if((length(gof_stat) == 1) && (gof_stat == "all")){
    gof_stat <- c("Geodesics","3-Motifs", "Dyad Shared Partners", "Edge Shared Partners", "Indegree","Outdegree","Degree","Incoming K-stars")
  }
  if(x$forms$directed & ("Degree"%in%gof_stat) & !(x$bipartite)){
    gof_stat <- c(gof_stat[-which(gof_stat=="Degree")],"Indegree","Outdegree")
  }
  if(x$bipartite & ("Degree"%in%gof_stat)){
    gof_stat <- c(gof_stat,"Degree Family 1", "Degree Family 2")
  }
  if(any(c("Outdegree","Indegree","3-Motifs")%in%gof_stat) & !x$forms$directed){
    stop("Requested statistic not meaningful for undirected networks.")
  }
  if(any(c("Dyad Shared Partners", "Edge Shared Partners","Incoming K-stars")%in%gof_stat)){
    if (!(requireNamespace("ergm", quietly = TRUE) & requireNamespace("network", quietly = TRUE))) {
      stop("Packages \"ergm\" and \"network\" needed to compute requested statistics. Please install them.",
           call. = FALSE)
    } else {
      message("Resorting to third-party package \"ergm\"; expect substantial increase in computation time.\n")
    }
  }
  
  ## Define helper function
  gof_getter <- function(gof, nets, fm){
    mat_maker <- ifelse(fm$bipartite, igraph::as_incidence_matrix,igraph::as_adjacency_matrix)
    switch(gof,
           "Indegree" = sapply(nets,
                               function(y){
                                 igraph::degree_distribution(y,
                                                             mode = "in")[1:(nrow(fm$monadic.data[[1]])-2)]
                               }),
           "Outdegree" = sapply(nets,
                                function(y){
                                  igraph::degree_distribution(y,
                                                              mode = "out")[1:(nrow(fm$monadic.data[[1]])-2)]
                                }),
           "Degree" = sapply(nets,
                             function(y){
                               if(x$bipartite){
                                 max_x <- max(nrow(fm$monadic.data[[1]]), 
                                              nrow(fm$monadic.data[[1]]))
                               } else {
                                 max_x <- nrow(fm$monadic.data[[1]])
                               }
                               igraph::degree_distribution(y,
                                                           mode = "all")[1:(max_x-2)]
                             }),
           "Degree Family 1" = sapply(nets,
                                      function(y){
                                        types <- igraph::vertex_attr(y, "type")
                                        igraph::degree_distribution(y, v = (igraph::V(y))[types])[1:(max(fm$NodeIndex[,2])-1)]
                                      }),
           "Degree Family 2" = sapply(nets,
                                      function(y){
                                        types <- !(igraph::vertex_attr(y, "type"))
                                        igraph::degree_distribution(y, v =( igraph::V(y))[types])[1:(max(fm$NodeIndex[,1])-1)]
                                      }),
           "Geodesics" = lapply(nets,
                                function(y){
                                  prop.table(table(igraph::distances(y)))
                                }),
           "3-Motifs" = sapply(nets,
                               function(y){
                                 res <- igraph::triad_census(y)
                                 res/sum(res, na.rm=TRUE)
                               }),
           "Dyad Shared Partners" = sapply(nets,
                                           function(y){
                                             prop.table(getS3method("summary", "formula", envir=asNamespace("ergm"))(network::network(mat_maker(y, sparse=FALSE), directed=igraph::is_directed(y)) ~
                                                                                                                       dsp(0:(igraph::gorder(y) - 2))))
                                           }),
           "Edge Shared Partners" = sapply(nets,
                                           function(y){
                                             prop.table(getS3method("summary", "formula", envir=asNamespace("ergm"))(network::network(mat_maker(y, sparse=FALSE), directed=igraph::is_directed(y)) ~
                                                                                                                       esp(0:(igraph::gorder(y) - 2))))
                                           }),
           "Incoming K-stars" = sapply(nets,
                                       function(y){
                                         prop.table(getS3method("summary", "formula", envir=asNamespace("ergm"))(network::network(mat_maker(y, sparse=FALSE), directed=igraph::is_directed(y)) ~
                                                                                                                   istar(0:(igraph::gorder(y) - 1))))
                                       })
    )
  }
  
  
  # Get networks
  
  el <- NetMix::simulate.mmsbmB(x, samples, seed=seed,
                                new.data.dyad,
                                new.data.monad)
  if(!is.null(new.data.dyad)){
    if(is.null(x$forms$timeID)){
      tid <- "(tid)"
      new.data.dyad[,tid] <- 1
    } else {
      tid <- x$forms$t_id_d
    }
    var_names <- c(with(x$forms, c(senderID, receiverID)), tid)
    new_y <- new.data.dyad[, all.vars(x$forms$formula.dyad)[1]]
    obs_dyad <- new.data.dyad[, var_names]
    el_obs <- new.data.dyad[new_y==1, var_names]
  } else {
    obs_dyad <- x$dyadic.data[,c("(sid)","(rid)","(tid)")]
    el_obs <- x$dyadic.data[x$Y==1,c("(sid)","(rid)", "(tid)")]
  }
  
  
  ## Convert to igraph objects
  nets_sim <- lapply(el,
                     function(i){
                       x_full <- cbind(obs_dyad, i)
                       x_sub <- x_full[x_full[,4] == 1, c(1, 2, 3)]
                       lapply(seq.int(length(unique(x_full[,3]))),
                              function(y){
                                x_sub_y <- x_sub[x_sub[,3]==y, c(1,2)]
                                if(x$bipartite){
                                  tmp_g <- igraph::make_empty_graph(directed = FALSE)
                                  tmp_g <- igraph::add_vertices(tmp_g,
                                                                nv = length(unique(x_sub_y[,1])),
                                                                attr = list(name = unique(x_sub_y[,1]),
                                                                            type = rep(0, length(unique(x_sub_y[,1])))))
                                  tmp_g <- igraph::add_vertices(tmp_g,
                                                                nv = length(unique(x_sub_y[,2])),
                                                                attr = list(name = unique(x_sub_y[,2]),
                                                                            type = 1))
                                  tmp_g <- igraph::add_edges(tmp_g, as.vector(t(x_sub_y)))
                                  return(tmp_g)
                                } else {
                                  return(igraph::graph_from_edgelist(as.matrix(x_sub_y), x$forms$directed))
                                }
                              })
                     })
  el_obs_list <- split.data.frame(el_obs, el_obs[,3])
  net_obs <- lapply(el_obs_list,
                    function(y){
                      if(x$bipartite){
                        tmp_g <- igraph::make_empty_graph(directed = FALSE)
                        tmp_g <- igraph::add_vertices(tmp_g,
                                                      nv = length(unique(y[,1])),
                                                      attr = list(name = unique(y[,1]),
                                                                  type = rep(0, length(unique(y[,1])))))
                        tmp_g <- igraph::add_vertices(tmp_g,
                                                      nv = length(unique(y[,2])),
                                                      attr = list(name = unique(y[,2]),
                                                                  type = 1))
                        tmp_g <- igraph::add_edges(tmp_g, as.vector(t(y[,c(1,2)])))
                        return(tmp_g)
                      } else {
                        return(igraph::graph_from_edgelist(as.matrix(y[,c(1, 2)]), x$forms$directed))
                      }
                    })
  
  ## Compute for simulated nets
  sim_stats_l <- lapply(gof_stat, gof_getter, nets = unlist(nets_sim, recursive = FALSE), fm = x)
  alpha <- (1 - level)/2
  sim_stats <- mapply(function(z, y){
    if(is.list(z)){
      z <- do.call(NetMix:::.cbind.fill, z)
    }
    z[is.na(z)] <- 0
    res <- as.data.frame(t(apply(z, 1, quantile, probs = c(alpha, 0.5, level + alpha), na.rm=TRUE)))
    names(res) <- c("LB","Est","UB")
    res$GOF <- y
    res$Val <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rownames(res)))
    return(res)
  },
  sim_stats_l, gof_stat,
  SIMPLIFY = FALSE)
  sim_stats_full <-  do.call("rbind",sim_stats)
  Observed_l <- lapply(gof_stat, gof_getter, nets = net_obs, fm = x)
  obs_stats <- mapply(function(z, y){
    if(is.list(z)){
      z <- do.call(NetMix:::.cbind.fill, z)
    }
    z[is.na(z)] <- 0
    res <- data.frame(Observed = apply(z, 1, median, na.rm=TRUE))
    res$GOF <- y
    res$Val <- as.numeric(gsub(".*?([0-9]+).*", "\\1", rownames(res)))
    return(res)
  },
  Observed_l, gof_stat,
  SIMPLIFY = FALSE)
  Observed <- do.call("rbind", obs_stats)
  res_df <- merge(sim_stats_full, Observed)
  
  ## Plot results
  return(ggplot2::ggplot(data=res_df, ggplot2::aes_string(x="Val", y="Observed")) +
           ggplot2::facet_wrap(~GOF, scales="free") +
           ggplot2::geom_linerange(ggplot2::aes_string(ymin="LB", ymax="UB"), col="black", lwd=2) +
           ggplot2::geom_line(ggplot2::aes_string(y="Observed"), lwd=1.1, alpha=0.8, col="darkred") +
           ggplot2::theme_bw() +
           ggplot2::xlab("") +
           ggplot2::ylab("Density"))
} 
