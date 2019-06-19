#' Plot posterior predictive checks using structural network charactericts
#'
#' The function generates a variety of plots that serve as posterior predictive checks on the goodness of fit of a fitted \code{mmsbm} object.
#'
#' Goodness of fit of network models has typically been established by evaluating how the structural characteristics of predicted networks 
#' compare to those of the observed network. When estimated in a Bayesian framework, this approach is equivalent to 
#' conducting posterior preditive checks on these structural quantities of interest. When \code{new.data.dyad} and/or \code{new.data.monad} are
#' passed that are different from those used in estimation, this is equivalent to conducting posterior predictive checks out-of-sample. 
#' 
#' The set of structural features used to determine goodness of fit is somewhat arbitrary, and chosen mostly to incorporate various 
#' first order, second order, and (to the extent possible) third-order characteristics of the network. "Geodesics" focuses on the distribution over 
#' observed and predicted geodesic distances between nodes; "Indegree" and "Outdegree" focuses on the distribution over incoming and outgoing connections 
#' per node; "3-motifs" focus on a distribution over possible connectivity patterns between triads; "Dyad Shared Partners" focuses on the distribution
#' over the number of shared partners between any two dayds; "Edge Shared Partners" is similarly defined, but w.r.t. edges, rather than dyads; and finally
#' "Incoming K-stars" focuses on a frequency distribution over stars with k=1,... spokes. 
#' 
#' Obtaining samples of the last three structural features can be very computationally expensive, and is discouraged on networks with more than 50
#'  nodes.
#'  
#'  
#' @param fm An object of class \code{mmsbm}, a result of a call to \code{mmsbm}.
#' @param gof_stat Character vector. Accepts any subset from "Geodesics", "3-Motifs", "Indegree", "Outdegree",
#'                 "Dyad Shared Partners", "Edge Shared Partners", and "Incoming K-stars". See details.
#' @param level Double. Level of credible interval for posterior predictive distribution around structural quantities of interest. 
#' @param samples Integer. Number of sampled networks from model's posterior predictive using \code{\link{simulate.mmsbm}}.
#' @param new.data.dyad See \code{\link{simulate.mmsbm}}.
#' @param new.data.monad See \code{\link{simulate.mmsbm}}.
#' @param parametric_mm See \code{\link{simulate.mmsbm}}.
#'

gof.mmsbm <- function(fm,
                      gof_stat = c("Geodesics","Degree"),
                      level = 0.95,
                      samples = 50,
                      new.data.dyad = NULL,
                      new.data.monad  = NULL, 
                      parametric_mm = FALSE
                      ){
  ## Define helper function
  gof_getter <- function(gof, nets, fm){
    switch(gof,
           "Indegree" = sapply(nets,
                             function(y){
                               igraph::degree_distribution(y,
                                                           mode = "in")[1:(nrow(fm$monadic.data)-2)]
                             }),
           "Outdegree" = sapply(nets,
                             function(y){
                               igraph::degree_distribution(y,
                                                           mode = "out")[1:(nrow(fm$monadic.data)-2)]
                             }),
           "Degree" = sapply(nets,
                                function(y){
                                  igraph::degree_distribution(y,
                                                              mode = "all")[1:(nrow(fm$monadic.data)-2)]
                                }),
           "Geodesics" = sapply(nets,
                               function(y){
                                 prop.table(igraph::distance_table(y, fm$directed)$res)
                               }),
           "3-Motifs" = sapply(nets,
                            function(y){
                              res <- igraph::triad_census(y)
                              res/sum(res, na.rm=TRUE)
                            }),
           "Dyad Shared Partners" = sapply(nets,
                          function(y){
                            prop.table(getS3method("summary", "formula", envir=asNamespace("ergm"))(network::network(igraph::as_adj(y, sparse=FALSE), directed=igraph::is_directed(y)) ~
                                         dsp(0:(igraph::gorder(y) - 2))))
                            }),
           "Edge Shared Partners" = sapply(nets,
                          function(y){
                            prop.table(getS3method("summary", "formula", envir=asNamespace("ergm"))(network::network(igraph::as_adj(y, sparse=FALSE), directed=igraph::is_directed(y)) ~
                                         esp(0:(igraph::gorder(y) - 2))))
                          }),
           "Incoming K-stars" = sapply(nets,
                            function(y){
                              prop.table(getS3method("summary", "formula", envir=asNamespace("ergm"))(network::network(igraph::as_adj(y, sparse=FALSE), directed=igraph::is_directed(y)) ~
                                           istar(0:(igraph::gorder(y) - 1))))
                            })
    )
  }

  
  ## Get networks
  el <- replicate(samples, simulate(fm, new.data.dyad,
                                    new.data.monad, 
                                    parametric_mm), simplify = FALSE) 
  if(!is.null(new.data.dyad)){
    if(is.null(fm$forms$timeID)){
      tid <- "(tid)"
      new.data.dyad[,tid] <- 1
    } else {
      tid <- fm$forms$tid
    }
    var_names <- c(with(fm$forms, c(senderID, receiverID)), tid)
    new_y <- new.data.dyad[, all.vars(net3.dmodel$forms$formula.dyad)[1]]
    obs_dyad <- new.data.dyad[, var_names]
    el_obs <- new.data.dyad[new_y==1, var_names]
  } else {
    obs_dyad <- fm$dyadic.data[,c("(sid)","(rid)","(tid)")]
    el_obs <- fm$dyadic.data[fm$Y==1,c("(sid)","(rid)", "(tid)")]
  }
  ## Convert to igraph objects
  nets <- lapply(el,
                 function(x){
                   x_full <- cbind(obs_dyad, x)
                   x_sub <- x_full[x_full[,4] == 1, c(1, 2, 3)]
                   lapply(seq.int(length(unique(x_full[,3]))),
                          function(y){
                            x_sub_y <- x_sub[x_sub[,3]==y, c(1,2)]
                            igraph::graph_from_edgelist(as.matrix(x_sub_y), fm$directed)
                          })
                 })
  el_obs_list <- split.data.frame(el_obs, el_obs[,3])
  net_obs <- lapply(el_obs_list,
                    function(y){
                      igraph::graph_from_edgelist(as.matrix(y[,c(1, 2)]), fm$directed)
                    })
  
  if((length(gof_stat) == 1) && (gof_stat == "all")){
    gof_stat <- c("Geodesics","3-Motifs", "Dyad Shared Partners", "Edge Shared Partners", "Indegree","Outdegree","Degree","Incoming K-stars")
  }
  if(fm$directed & ("Degree"%in%gof_stat)){
    gof_stat <- c(gof_stat,"Indegree","Outdegree")
  }
  if(any(c("Outdegree","Indegree","3-Motifs")%in%gof_stat) & !fm$directed){
    stop("Requested statistic not meaningful for undirected networks.")
  }
  if(any(c("Dyad Shared Partners", "Edge Shared Partners","Incoming K-stars")%in%gof_stat)){
    cat("Resorting to third-party package `ergm'; expect substantial increase in computation time.\n")
  }
  
  ## Compute for simulated nets
  sim_stats_l <- lapply(gof_stat, gof_getter, nets = unlist(nets, recursive = FALSE), fm = fm)
  alpha <- (1 - level)/2 
  sim_stats <- mapply(function(x, y){
                       if(is.list(x)){
                         require(rowr)
                         x <- do.call(cbind.fill, x)
                       }
                       x[is.na(x)] <- 0  
                       res <- as.data.frame(t(apply(x, 1, quantile, probs = c(alpha, 0.5, level + alpha), na.rm=TRUE)))
                       names(res) <- c("LB","Est","UB")
                       res$GOF <- y
                       res$Val <- 1:nrow(res)
                       return(res)
                      },
                      sim_stats_l, gof_stat,
                      SIMPLIFY = FALSE)
  sim_stats_full <-  do.call("rbind",sim_stats)               
  Observed <- apply(do.call(rbind,lapply(gof_stat, gof_getter, nets = net_obs, fm = fm)),1,median)
  Observed[is.na(Observed)] <- 0
  sim_stats_full$Observed <- Observed

  ## Plot results
  ggplot2::ggplot(subset(sim_stats_full, Est > 0), aes(x=Val, y=Observed)) +
    ggplot2::facet_wrap(~GOF, scales="free") +
    ggplot2::geom_linerange(aes(ymin=LB, ymax=UB), col="gray60", lwd=2) +
    ggplot2::geom_line(aes(y=Observed), lwd=1.1) +
    ggplot2::theme_bw() + 
    ggplot2::xlab("") +
    ggplot2::ylab("Density")
} 