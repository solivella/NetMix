#########################
## Misc. Helper functions
#########################

gof <- function (x, ...) {
  UseMethod("gof", x)
}

.cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(NA, n-nrow(x), ncol(x))))) 
}

.mpower <- function(x, p){
  orig <- x
  while(p > 1){
    x <- x %*% orig
    p = p - 1 
  }
  return(x)
}

.reRange <- function(x) {
  nmax <- 1.0 - 1e-12
  nmin <- 1e-12
  (nmax - nmin) * (x - 1.0) + nmax
}
.modSum <- function(x, y){
  replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
}

.bernKL <- function(x,y){
  x * log(x/y) + (1-x)*log((1-x)/(1-y))
}

## Permute blockmodels
.findPerm <- function(block.list, target.mat=NULL, use.perms=TRUE){
  tar_nprov <- is.null(target.mat) 
  if(tar_nprov){
    target.mat <- block.list[[1]] 
  }
  n <- ncol(target.mat)
  if(use.perms){
    require(Matrix, quietly = TRUE); require(gtools, quietly = TRUE) 
    all_perms <- gtools::permutations(n, n)
  } else {
    require(igraph, quietly = TRUE);
  }
  res <- lapply(seq_along(block.list),
                function(x){
                  if(use.perms){
                    norms <- apply(all_perms, 1,
                                   function(p){
                                     base::norm(block.list[[x]][p, p]-target.mat, type="f")
                                   })
                    P <- as.matrix(as(all_perms[which.min(norms),], "pMatrix"))
                  } else {
                    P <- as.matrix(igraph::match_vertices(plogis(block.list[[x]]),
                                                          plogis(target.mat),
                                                          m = 0,
                                                          start = diag(ncol(block.list[[x]])),
                                                          iteration = 10)$P)
                  }
                  if(tar_nprov){
                    target.mat <<- t(P) %*% block.list[[x]] %*% P
                  }
                  return(P)
                })
  return(res)
}

.findPerm2 <- function(block.list, target.mat=NULL, use.perms=TRUE){
  nblocks <- ncol(block.list[[1]])
  if(use.perms){
    require(Matrix, quietly = TRUE); require(gtools, quietly = TRUE)
    all_perms <- gtools::permutations(nblocks, nblocks)
    n.perms <- nrow(all_perms)
    internal_loss <- vector("list", n.perms)
  } else {
    require(igraph, quietly = TRUE)
  }
  nsteps <- length(block.list)
  for(i in 1:nsteps){
    block.list[[i]] <- (block.list[[i]])
  }
  if(is.null(target.mat)){
    target.mat <- Reduce("+", block.list)/length(block.list)
    nulltarg <- TRUE 
  } else {
    nulltarg <- FALSE
  }
  new_loss <- sum(sapply(block.list, 
                         function(mat, target = target.mat){
                           base::norm(mat-target, type="o")
                           #sum(block.list[[i]] * log(block.list[[i]]/target.mat))
                         }))
  delta <- 1e6
  iter <- 1
  perms <- vector("list", nsteps)
  losses <- vector("list", nsteps)
  while((iter < 100) & (delta > 1e-2)){
    old_loss <- new_loss
    ## step 1
    if(nulltarg){
      target.mat <- Reduce("+", block.list)/length(block.list)
    }
    ## step 2
    old.list <- block.list
    for(i in 1:nsteps){
      if(use.perms){
        norms <- apply(all_perms, 1,
                       function(p){
                         base::norm(block.list[[i]][p, p]-target.mat, type="o")
                         #sum(block.list[[i]][p, p] * log(block.list[[i]][p, p]/target.mat))
                       })
        perms[[i]] <- as.matrix(as(all_perms[which.min(norms),], "pMatrix"))
      } else {
        perms[[i]] <- as.matrix(igraph::match_vertices(block.list[[i]], target.mat, m = 0, start = diag(nblocks), 100)$P)
      }
      block.list[[i]] <- t(perms[[i]]) %*% block.list[[i]] %*% perms[[i]]
      losses[[i]] <- base::norm(block.list[[i]] - target.mat, type="o")
      #losses[[i]] <- sum(block.list[[i]] * log(block.list[[i]]/target.mat))
    }
    ## compute loss and delta
    new_loss <- Reduce("+", losses)
    delta <- abs((new_loss - old_loss)/old_loss)
    print(delta)
    iter <- iter + 1
  }
  return(perms)
}




.transf <- function(mat){
  (mat * (nrow(mat) - 1) + 1/ncol(mat))/nrow(mat)
}
.mapID <- function(uvec, vec){
  temp_map <- seq_along(uvec)
  names(temp_map) <- as.character(uvec)
  temp_map[as.character(vec)]
}


cluster.mems <- function(fm, t=unique(fm$dyadic.data$`(tid)`), n=10, demean=FALSE){
  Mem <- fm$MixedMembership[,fm$monadic.data[,"(tid)"] %in% t]
  Nodes <- unlist(lapply(strsplit(colnames(Mem), "@"), "[[", 1))
  node.mems <- t(do.call(cbind, lapply(unique(Nodes), function(x){
    rowMeans(as.matrix(Mem[,Nodes==x]))})))
  rownames(node.mems) <- as.character(unique(Nodes))
  if(demean){
    node.mems2 <- apply(node.mems, 2, function(x){x - mean(x)})
    sort(apply(node.mems2, 1, which.max))
  } else {
    lapply(1:fm$n_blocks, function(x){
      node.mems[order(node.mems[,x], decreasing=T)[1:n],x]})
  }
}




cluster.time.node <- function(fm){
  require(ggplot2)
  for(i in unique(fm$monadic.data[,"(nid)"])){
    avgmem <- lapply(1:nrow(fm$MixedMembership), function(x){fm$MixedMembership[x,which(fm$monadic.data[,"(nid)"]==i)]})
    pds <- as.character(fm$monadic.data[fm$monadic.data[,"(nid)"]==i, "(tid)"])
    avgmem <- as.data.frame(cbind(rep(pds, nrow(fm$MixedMembership)),
                                  unlist(avgmem),
                                  rep(1:nrow(fm$MixedMembership), each=length(pds))))
    colnames(avgmem) <- c("Time", "Membership", "Group")
    avgmem$Group <- factor(avgmem$Group, levels=4:1)
    if(class(avgmem$Membership) == "factor"){avgmem$Membership <- as.numeric(as.character(avgmem$Membership))}
    if(class(avgmem$Time) == "factor"){avgmem$Time <- as.numeric(as.character(avgmem$Time))}
    print(ggplot() + ggtitle(paste("Group Membership Over Time, Node", i)) + theme(plot.title = element_text(hjust = 0.5)) +
            geom_area(aes(y = Membership, x = Time, fill=Group), data = avgmem,
                      stat="identity", position="stack")  + guides(fill=guide_legend(title="Group")))
  }
}




cluster.time <- function(fm){
  avgmems <- lapply(1:nrow(fm$MixedMembership), function(x){
    tapply(fm$MixedMembership[x,], fm$monadic.data[,"(tid)"], mean)})
  avgmems <- as.data.frame(cbind(rep(unique(as.character(fm$monadic.data[,"(tid)"])), nrow(fm$MixedMembership)),unlist(avgmems),
                                 rep(1:nrow(fm$MixedMembership), each=length(unique(fm$monadic.data[,"(tid)"])))))
  colnames(avgmems) <- c("Time", "Avg.Membership", "Group")
  avgmems$Group <- factor(avgmems$Group, levels=length(unique(avgmems$Group)):1)
  if(class(avgmems$Avg.Membership) == "factor"){avgmems$Avg.Membership <- as.numeric(as.character(avgmems$Avg.Membership))}
  if(class(avgmems$Time) == "factor"){avgmems$Time <- as.numeric(as.character(avgmems$Time))}
  require(ggplot2)
  ggplot() + ggtitle("Average Group Membership Over Time") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_area(aes(y = Avg.Membership, x = Time, fill=Group), data = avgmems,
              stat="identity", position="stack")  + guides(fill=guide_legend(title="Group"))
}




cluster.node <- function(fm){
  node.mems <- t(do.call(cbind, lapply(unique(fm$monadic.data[,"(nid)"]), function(x){
    rowMeans(as.matrix(fm$MixedMembership[,fm$monadic.data[,"(nid)"]==x]))})))
  rownames(node.mems) <- as.character(unique(fm$monadic.data[,"(nid)"]))
  node.mems <- node.mems[order(node.mems[,1], decreasing=T),]
  require(plotly, quietly=TRUE)
  plot_ly(z=node.mems, type="heatmap", x=paste("Group", 1:ncol(node.mems)),
          y=rownames(node.mems), colorscale="Portland")
}



est.edgeplot <- function(fm){
  estimated <- predict(fm)
  obs.node <- sapply(unique(fm$dyadic.data[,"(sid)"]), function(x){
    sum(fm$dyadic.data[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x,1])})
  est.node <- sapply(unique(fm$dyadic.data[,"(sid)"]), function(x){
    sum(estimated[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x])})
  
  max.edges <- round(max(c(est.node, obs.node)))
  plot(obs.node, est.node, main = "Observed vs. Estimated Number of Edges per Node",
       xlim = c(round(min(c(est.node, obs.node),na.rm=T)), round(max(c(est.node, obs.node),na.rm=T))),
       ylim=c(round(min(c(est.node, obs.node),na.rm=T)), round(max(c(est.node, obs.node),na.rm=T))),
       xlab = "Observed Number of Edges",
       ylab = "Estimated Number of Edges")
}



degree.dist <- function(fm, Y){
  estimated <- predict(fm)
  obs.node <- sapply(unique(fm$dyadic.data[,"(sid)"]), function(x){
    sum(fm$dyadic.data[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x,1])})
  est.node <- sapply(unique(fm$dyadic.data[,"(sid)"]), function(x){
    sum(estimated[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x])})
  dens.o <- density(obs.node)
  dens.e <- density(est.node)
  plot(NA, xlab = "", ylab="", main="Degree Distribution",
       xlim=range(c(dens.o$x, dens.e$x)), ylim=range(c(dens.o$y, dens.e$y)))
  lines(dens.o)
  lines(dens.e, col="red")
  legend("topright", legend=c("Observed", "Estimated"),
         col=c("black", "red"), lty=1)
}


.pi.hat <- function(X, beta){
  if(dim(beta)[3]==1){
    mu <- exp(X %*% beta[,,1])
    pi.states <- list(t(mu))
  } else {
    pi.states <- lapply(1:dim(beta)[3], function(m){
      mu <- exp(X %*% beta[,,m])
      return(t(mu))
    })
  }
  return(pi.states)
}

.e.pi <- function(pi_l, kappa){
  if(is.null(dim(kappa))){
    return(pi_l[[1]])
  } else {
    n_states <- nrow(kappa)
    pi.states <- lapply(1:nrow(kappa),
                        function(m){
                          pi_l[[m]] * kappa[m,]
                        })
    return(Reduce("+", pi_l))
  }
}

## covFX <- function(fm, cov, shift, max.val=FALSE){
##   predict.ties <- predict(fm)
##   monadic.data2 <- fm$monadic.data
##   monadic.data2[,cov] <- fm$monadic.data[,cov] + shift
##   if(!isFALSE(max.val)){
##     monadic.data2[which(fm$monadic.data[,cov] == max(fm$monadic.data[,cov])),cov] <- max.val
##   }
##   predict.ties2 <- predict(fm, monad=monadic.data2)
##   FX <- list(mean(predict.ties2 - predict.ties), #avg
##              tapply(predict.ties2-predict.ties, fm$dyadic.data[,"(tid)"], mean), #time
##              sapply(unique(fm$monadic.data[,"(nid)"]), function(x){ #node
##                mean((predict.ties2-predict.ties)[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x])}),
##              tapply(predict.ties2-predict.ties, paste(fm$dyadic.data[,"(sid)"], fm$dyadic.data[,"(rid)"], sep="_"), mean),#dyad
##              predict.ties2 - predict.ties) #dyad-time
##   names(FX[[5]]) <- paste(fm$dyadic.data[,"(sid)"], fm$dyadic.data[,"(rid)"], sep="_")
##   names(FX) <- c(paste("Overall Avg. Effect of", cov), paste("Avg. Effect of", cov, "by Time"),
##                  paste("Avg. Effect of", cov, "by Node"), paste("Avg. Effect of", cov, "by Dyad"),
##                  paste("Effect of", cov, "by Dyad-Time"))
##   return(FX)
## }



boot.mmsbm <- function(fm, Y, cov, m=1000){  ## add multiple core function here (for each loop)
  results <- list()
  results.cov <- list()
  if(ncol(cov) != ncol(fm$coef)-1){cov <- cov[,match(colnames(fm$coef)[2:ncol(fm$coef)], colnames(cov))]}
  orig.order <- colSums(fm$pi) # max, min, middle (will need to adjust this)
  for(i in 1:m){
    samp <- sample(1:nrow(Y), nrow(Y), replace=TRUE)
    Y.samp <- Y[samp, samp]
    cov.samp <- cov[samp,]
    ## Adjust fit to re-do exact model, including formula, control list, and priors
    fit <- hmmsb(Y.samp, formula = ~., data = cov.samp,
                 n.blocks = ncol(fm$pi),
                 a = a_prior, b = b_prior, self.tie = 1,
                 cvb0=TRUE, mmsb.control = list(vi.iter=2000, em.iter=20, vi.tol = 1e-4))
    first <- which(colSums(fit$pi)==max(colSums(fit$pi)))
    second <- which(colSums(fit$pi)==min(colSums(fit$pi)))
    third <- (1:3)[-which(1:3 %in% c(first, second))]
    results[[i]] <- fit$coef[c(first, second, third),]
    #results.cov[[i]] <- covFX(fit, cov.samp, 1)[[1]]
  }
  return(results, results.cov)
}

## modularized fns from mmsbm.R

#Missing data
.missing <- function(formula.dyad, formula.monad1, formula.monad2=NULL, 
                     data.dyad, data.monad1, data.monad2=NULL, missing){
  
  if(missing=="indicator method"){
    # dyadic dataset
    if(length(all.vars(formula.dyad[[3]]))){
      miss.d <- apply(data.dyad[,all.vars(formula.dyad[[3]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.dyad)
      md <- names(miss.d[miss.d])
      if(length(md)>0){
        m.ind <- apply(as.data.frame(data.dyad[,md]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(md, "_missing", sep="")
        data.dyad[,md] <- apply(as.data.frame(data.dyad[,md]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        })
        data.dyad <- cbind(data.dyad, m.ind)
        fc <- paste(as.character(formula.dyad[[2]]), as.character(formula.dyad[[1]]),
                    paste(c(all.vars(formula.dyad)[-1], colnames(m.ind)), collapse=" + "))
        formula.dyad <- eval(parse(text=fc))
      }
    }
    # monadic1 dataset
    if(length(all.vars(formula.monad1[[2]]))){
      miss.m <- apply(data.monad1[,all.vars(formula.monad1[[2]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.monad1)
      mm <- names(miss.m[miss.m])
      if(length(mm)>0){
        m.ind <- apply(as.data.frame(data.monad1[,mm]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(mm, "_missing", sep="")
        data.monad1[,mm] <- as.vector(apply(as.data.frame(data.monad1[,mm]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        }))
        data.monad1 <- cbind(data.monad1, m.ind)
        fc <- paste("~", paste(c(all.vars(formula.monad1), colnames(m.ind)),  collapse=" + "))
        formula.monad1 <- eval(parse(text=fc))
      }
    }
    # monadic2 dataset
    if(length(all.vars(formula.monad2[[2]]))){
      miss.m <- apply(data.monad2[,all.vars(formula.monad2[[2]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.monad2)
      mm <- names(miss.m[miss.m])
      if(length(mm)>0){
        m.ind <- apply(as.data.frame(data.monad2[,mm]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(mm, "_missing", sep="")
        data.monad2[,mm] <- as.vector(apply(as.data.frame(data.monad2[,mm]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        }))
        data.monad2 <- cbind(data.monad2, m.ind)
        fc <- paste("~", paste(c(all.vars(formula.monad2), colnames(m.ind)),  collapse=" + "))
        formula.monad2 <- eval(parse(text=fc))
      }
    }
  }
  if(missing=="listwise deletion"){
    if(length(all.vars(formula.dyad[[3]]))){
      mdyad <- apply(data.dyad[,all.vars(formula.dyad[[3]])], 1, function(x){!any(is.na(x))})
    } else {
      mdyad <- TRUE
    }
    #monad1
    if(length(all.vars(formula.monad1[[2]]))){
      mmonad1 <- apply(data.monad1[,all.vars(formula.monad1[[2]])], 1, function(x){!any(is.na(x))})
    } else {
      mmonad1 <- TRUE
    }
    #monad2
    if(length(all.vars(formula.monad2[[2]]))){
      mmonad2 <- apply(data.monad2[,all.vars(formula.monad2[[2]])], 1, function(x){!any(is.na(x))})
    } else {
      mmonad2 <- TRUE
    }
    
    data.dyad <- data.dyad[mdyad,]
    data.monad1 <- data.monad1[mmonad1,]
    data.monad2 <- data.monad2[mmonad2,]
    d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
      nts <- data.monad1[data.monad1[,timeID]==x,nodeID] #### How to change this to also include data.monad2?
      dd <- data.dyad[data.dyad[,timeID]==x,]
      dd <- dd[dd[,senderID] %in% nts & dd[,receiverID] %in% nts,]
      return(dd)
    })
    data.dyad <- do.call("rbind", d.keep)
  }
}

#Create sociomatrix
.createSocio <- function(dyads, all.nodes1, all.nodes2){
res <- lapply(dyads,
                   function(dyad_df,
                            nnode1 = length(all.nodes1),
                            nnode2 = length(all.nodes2),
                            nodes1 = all.nodes1,
                            nodes2 = all.nodes2)
#if bipartite:
#need adj_mat to be composed of nnode1=length(all.nodes1), nnode2=length(all.nodes2), and dimnames are nodes1, nodes2
                     # then create two sets of node_names
                     #create function here that takes 1 and 2 always; then double up in input ifit is bipartite (outside this fn, so remove bipartite arg)
                   {
                     indeces <- as.matrix(dyad_df[,c("(sid)","(rid)")])
                     time_ind <- unique(dyad_df[,"(tid)"])
                     adj_mat <-  matrix(NA,
                                        nnode1,
                                        nnode2,
                                        dimnames = list(nodes1,
                                                        nodes2))
                     adj_mat[indeces] <- dyad_df[,4] # out of bounds
                     if(!directed){
                       adj_mat[indeces[,c(2,1)]] <- dyad_df[,4]
                     }
                     diag(adj_mat) <- 0
                     adj_mat[is.na(adj_mat)] <- sample(0:1, sum(is.na(adj_mat)), replace = TRUE)
                     if(!directed){
                       mat_ind <- which(upper.tri(adj_mat), arr.ind = TRUE)
                       adj_mat[mat_ind[,c(2,1)]] <- adj_mat[upper.tri(adj_mat)]
                     }
                     node_names1 <- paste(nodes1, "@", time_ind, sep="")
                     node_names2 <- paste(nodes2, "@", time_ind, sep="")
                     dimnames(adj_mat) <- list(node_names1,
                                               node_names2)
                     return(adj_mat)
                   })
return(res)
}

#' @rdname auxfuns
.initPi <- function(soc_mats,
                    Y,
                    dyads,
                    edges,
                    t_id_d,
                    t_id_n1,t_id_n2,
                    nodes_pp,nodes_pp1,nodes_pp2,
                    dyads_pp,
                    nt_id1,nt_id2,
                    node_id_period1,node_id_period2,
                    mu_b,
                    var_b,
                    n_dyads, n.blocks1, n.blocks2, periods, directed, ctrl){
  max_ll <- -Inf
  best_pi1 <- matrix(as.double(0.0), ncol=n_dyads, nrow=n.blocks1)
  best_pi2 <- matrix(as.double(0.0), ncol=n_dyads, nrow=n.blocks2)
  rownames(best_pi1) <- 1:n.blocks1
  rownames(best_pi2) <- 1:n.blocks2
  for(n in 1:ctrl$nstarts){
    temp_res <- vector("list", periods)
    for(i in 1:periods){
      if(!ctrl$init_gibbs) {
        mn <- ncol(soc_mats[[i]]) 
        if(directed){
          D_o <- 1/sqrt(.rowSums(soc_mats[[i]], mn, mn) + 1)
          D_i <- 1/sqrt(.colSums(soc_mats[[i]], mn, mn) + 1)
          C_o <- t(D_o * soc_mats[[i]])
          C_i <- t(D_i * t(soc_mats[[i]]))
          U <- t(C_o * D_i) %*% C_o +
            t(C_i * D_o) %*% C_i
        } else {
          D <- 1/sqrt(.rowSums(soc_mats[[i]], mn, mn) + 1)
          U <- t(D * soc_mats[[i]]) * D
        }
        if(ctrl$spectral) {
          n_elem1 <- n.blocks1 + 1
          res1 <- RSpectra::eigs_sym(U, n_elem1)
          eta1 <- res1$vectors[,1:n_elem1] %*% diag(res1$values[1:n_elem1])
          target1 <- eta1[,2:n_elem1] / (eta1[,1] + 1e-8)
          sig1 <- 1 - res1$values[n_elem1] / (res1$values[n.blocks1])
          sig1 <- ifelse(is.finite(sig1), sig1, 0)
          if(sig1 > 0.1){
            target1 <- target1[,1:(n_elem1 - 2), drop = FALSE]
          }
          n_elem2 <- n.blocks2 + 1
          res2 <- RSpectra::eigs_sym(t(U), n_elem2) ##transpose U
          eta2 <- res2$vectors[,1:n_elem2] %*% diag(res2$values[1:n_elem2])
          target2 <- eta2[,2:n_elem2] / (eta2[,1] + 1e-8)
          sig2 <- 1 - res2$values[n_elem2] / (res2$values[n.blocks2])
          sig2 <- ifelse(is.finite(sig2), sig2, 0)
          if(sig2 > 0.1){
            target2 <- target2[,1:(n_elem2 - 2), drop = FALSE]
          }
        } else {
          target1 <- U
          target2 <- t(U)
        }
        if(nrow(unique(target1)) > n.blocks1){
          clust_internal1 <- fitted(kmeans(x = target1,
                                          centers = n.blocks1,
                                          iter.max = 15,
                                          nstart = 1), "classes")
          
        } else {
          init_c1 <- sample(1:nrow(target1), n.blocks1, replace = FALSE)
          cents1 <- jitter(target[init_c1, ])
          clust_internal1 <- fitted(suppressWarnings(kmeans(x = target1,
                                                           centers = cents1,
                                                           iter.max = 15,
                                                           algorithm = "Lloyd",
                                                           nstart = 1)), "classes")
        }
        if(nrow(unique(target2)) > n.blocks2){
          clust_internal2 <- fitted(kmeans(x = target2,
                                           centers = n.blocks2,
                                           iter.max = 15,
                                           nstart = 1), "classes")
          
        } else {
          init_c2 <- sample(1:nrow(target2), n.blocks2, replace = FALSE)
          cents2 <- jitter(target[init_c2, ])
          clust_internal2 <- fitted(suppressWarnings(kmeans(x = target2,
                                                            centers = cents2,
                                                            iter.max = 15,
                                                            algorithm = "Lloyd",
                                                            nstart = 1)), "classes")
        }
        
        pi_internal1 <- model.matrix(~ factor(clust_internal1, 1:n.blocks1) - 1)
        pi_internal2 <- model.matrix(~ factor(clust_internal2, 1:n.blocks2) - 1)
        pi_internal1 <- .transf(pi_internal1)
        pi_internal2 <- .transf(pi_internal2)
        rownames(pi_internal1) <- rownames(soc_mats[[i]])
        rownames(pi_internal2) <- colnames(soc_mats[[i]])
        colnames(pi_internal1) <- 1:n.blocks1
        colnames(pi_internal2) <- 1:n.blocks2
        MixedMembership1 <- t(pi_internal1)
        MixedMembership2 <- t(pi_internal2)
        #int_dyad_id <- apply(dyads[[i]], 2, function(x)match(x, colnames(MixedMembership)) - 1)
        int_dyad_id <- dyads[[i]][,names(dyads[[i]])%in%c("(sid)","(rid)")]
        BlockModel <- approxB(edges[[i]], as.matrix(int_dyad_id), MixedMembership1, MixedMembership2)
        temp_res[[i]] <- list(BlockModel = BlockModel,
                              MixedMembership1 = MixedMembership1, MixedMembership2 = MixedMembership2)                         
        
      } else {
        #MM1
        n_prior <- (dyads_pp[i] - nodes_pp1[i]) * .05
        a <- plogis(ctrl$mu_b) * n_prior
        b <- n_prior - a
        lda_beta_prior <- lapply(list(b,a),
                                 function(prior){
                                   mat <- matrix(prior[2], n.blocks1, n.blocks1)
                                   diag(mat) <- prior[1]
                                   return(mat)
                                 })
        ret <- lda::mmsb.collapsed.gibbs.sampler(network = soc_mats[[i]],
                                                 K = n.blocks1,
                                                 num.iterations = 100L,
                                                 burnin = 50L,
                                                 alpha = ctrl$alpha,
                                                 beta.prior = lda_beta_prior)
        MixedMembership1 <- prop.table(ret$document_expects, 2)
        #MM2
        n_prior <- (dyads_pp[i] - nodes_pp2[i]) * .05
        a <- plogis(ctrl$mu_b) * n_prior
        b <- n_prior - a
        lda_beta_prior <- lapply(list(b,a),
                                 function(prior){
                                   mat <- matrix(prior[2], n.blocks2, n.blocks2)
                                   diag(mat) <- prior[1]
                                   return(mat)
                                 })
        ret <- lda::mmsb.collapsed.gibbs.sampler(network = soc_mats[[i]],
                                                 K = n.blocks2,
                                                 num.iterations = 100L,
                                                 burnin = 50L,
                                                 alpha = ctrl$alpha,
                                                 beta.prior = lda_beta_prior)
        MixedMembership2 <- prop.table(ret$document_expects, 2)
        
        
        colnames(MixedMembership1) <- rownames(soc_mats[[i]])
        colnames(MixedMembership2) <- colnames(soc_mats[[i]])
        
        #int_dyad_id <- apply(dyads[[i]], 2, function(x)match(x, colnames(MixedMembership)) - 1)
        int_dyad_id <- dyads[[i]][,names(dyads[[i]])%in%c("(sid)","(rid)")]
        BlockModel <- approxB(edges[[i]], as.matrix(int_dyad_id), MixedMembership1, MixedMembership2)
        if(any(is.nan(BlockModel))){
          BlockModel[is.nan(BlockModel)] <- 0.0
        }
        temp_res[[i]] <- list(BlockModel = BlockModel,
                              MixedMembership1 = MixedMembership1, MixedMembership2 = MixedMembership2)
        
      }
    }
    block_models <- lapply(temp_res, function(x)x$BlockModel)
    target_ind1 <- which.max(sapply(soc_mats, nrow))
    target_ind2 <- which.max(sapply(soc_mats, ncol))
    
    perms_temp1 <- .findPerm(block_models, target_mat = block_models[[target_ind1]], use_perms = ctrl$permute)
    perms_temp2 <- .findPerm(block_models, target_mat = block_models[[target_ind2]], use_perms = ctrl$permute)
    pis_temp1 <- lapply(temp_res, function(x)x$MixedMembership1) 
    pis_temp2 <- lapply(temp_res, function(x)x$MixedMembership2) 
    pi.ord1 <- as.numeric(lapply(pis_temp1, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order ## CHECK
    pi.ord2 <- as.numeric(lapply(pis_temp2, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order
    
    pi_init_t1 <- do.call(cbind,mapply(function(phi_tmp1,perm1){perm1 %*% phi_tmp1},
                                      pis_temp1[order(pi.ord1)], perms_temp1, SIMPLIFY = FALSE)) 
    pi_init_t2 <- do.call(cbind,mapply(function(phi_tmp2,perm2){perm2 %*% phi_tmp2},
                                       pis_temp2[order(pi.ord2)], perms_temp2, SIMPLIFY = FALSE)) 
    
    rownames(pi_init_t1) <- 1:n.blocks1
    rownames(pi_init_t2) <- 1:n.blocks2
    Z_tmp <- matrix(0, ncol = n_dyads, nrow = 1)
    X1_tmp <- matrix(1, ncol = ncol(pi_init_t1), nrow = 1)
    X2_tmp <- matrix(1, ncol = ncol(pi_init_t2), nrow = 1)
    ctrl$vi_iter <- 2
    ctrl$verbose <- FALSE
    ##sparsity <- mean(Y >= 0.5)
    fit_tmp <- mmsbm_fit(Z_tmp,
                         ##Z_tmp,
                         X1_tmp,
                         X2_tmp,
                         Y,
                         ##Y,
                         t_id_d,
                         t_id_n1,t_id_n2,
                         nodes_pp,nodes_pp1,nodes_pp2,
                         ##nt_id,
                         nt_id1,nt_id2,
                         node_id_period1,node_id_period2,
                         mu_b,
                         var_b,
                         array(0.0, c(1, n.blocks1, ctrl$states)),
                         array(1.0, c(1, n.blocks1, ctrl$states)),
                         array(0.0, c(1, n.blocks1, ctrl$states)),
                         array(1.0, c(1, n.blocks1, ctrl$states)),
                         0.0,
                         1.0,
                         pi_init_t1,pi_init_t2,
                         ctrl$kappa_init_t,
                         qlogis(block_models[[target_ind1]]), ##currently only using 1 family
                         array(ctrl$alpha, c(1,n.blocks1,ctrl$states)),array(ctrl$alpha, c(1,n.blocks2,ctrl$states)),
                         ctrl$gamma_init,
                         0.0,
                         ##sparsity,
                         ctrl)
    if(fit_tmp$LowerBound >= max_ll){
      best_pi1 <- pi_init_t1
      best_pi2 <- pi_init_t2
    }
  }
  return(list(best_pi1,best_pi2))
}

#transparency of colors for plotwebB
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
