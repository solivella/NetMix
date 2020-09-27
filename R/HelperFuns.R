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

# #Handle null phi
# .nullPhi <- function(
#   # full data passing through
#                      sociomatrix, #soc_mats
#                      ntid,
#                      mfm,
#                      dyads_pp,
#                      nodes_pp,
#                      state_init,
#                      #other items
#                      blocks,
#                      family2,
#                      periods,
#                      #control items
#                      alpha,
#                      conv_tol,
#                      em_iter,
#                      init.dyn.gibbs,
#                      mu_b,
#                      permute,
#                      phi_init_orig, #ctrl$phi_init_t
#                      seed,
#                      spectral,
#                      threads,
#                      var_b
# )
# {
#   offsetval = 0.0
#   if(family2){offsetval=length(sociomatrix)}
#   
#   if(is.null(phi_init_orig) & (periods == 1)){
#     phi_init_temp <- lapply(sociomatrix, function(mat){
#       if(directed){
#         D_o <- 1/sqrt(rowSums(mat) + 1)
#         D_i <- 1/sqrt(colSums(mat) + 1)
#         C_o <- t(D_o * mat)
#         C_i <- t(D_i * t(mat))
#         U <- t(C_o * D_i) %*% C_o +
#           t(C_i * D_o) %*% C_i
#       } else {
#         D <- 1/sqrt(rowSums(mat) + 1)
#         U <- t(D * mat) * D
#       }
#       if(spectral) {
#         n_elem <- blocks + 1
#         res <- RSpectra::eigs_sym(U, n_elem)
#         eta <- res$vectors[,1:n_elem] %*% diag(res$values[1:n_elem])
#         target <- eta[,2:n_elem] / (eta[,1] + 1e-8)
#         sig <- 1 - res$values[n_elem] / (res$values[blocks])
#         sig <- ifelse(is.finite(sig), sig, 0)
#         if(sig > 0.1){
#           target <- target[,1:(n_elem - 2)]
#         }
#       } else {
#         target <- U
#       }
#       init_c <- sample(1:nrow(target), blocks, replace = FALSE)
#       clust_internal <- fitted(kmeans(target,
#                                       blocks,
#                                       centers = target[init_c,],
#                                       #algorithm = "Lloyd",
#                                       nstart = 10), "classes")
#       # clust_internal <- cluster::clara(target,
#       #                                  blocks,
#       #                                  samples = 10,
#       #                                  sampsize = min(nrow(target), 100 + 2 * blocks),
#       #                                  rngR = TRUE,
#       #                                  correct.d = TRUE)$clustering
#       
#       phi_internal <- model.matrix(~ as.factor(clust_internal) - 1)
#       phi_internal <- .transf(phi_internal)
#       rownames(phi_internal) <- rownames(mat)
#       colnames(phi_internal) <- 1:blocks
#       t(phi_internal)
#     })
#     #phi_init_orig <- do.call(cbind, phi_init_temp)
#     #ctrl$phi_list <- phi_init_temp
#     res <- do.call(cbind, phi_init_temp)
#   }
#   else if (is.null(phi_init_orig) & (periods > 1)){
#     temp_res <- vector("list", periods)
#     for(i in 1:periods){
#       if (init.dyn.gibbs) {
#         n_prior <- (dyads_pp[i] - nodes_pp[offsetval+i]) * .05
#         a <- plogis(mu_b) * n_prior
#         b <- n_prior - a
#         lda_beta_prior <- lapply(list(b,a),
#                                  function(prior){
#                                    mat <- matrix(prior[2], blocks, blocks) # family1: n.blk1/n.blk1
#                                    diag(mat) <- prior[1]
#                                    return(mat)
#                                  })
#         ret <- lda::mmsb.collapsed.gibbs.sampler(network = sociomatrix[[i]],
#                                                  K = blocks,
#                                                  num.iterations = 100L,
#                                                  burnin = 50L,
#                                                  alpha = alpha, #alpha can be 2x
#                                                  beta.prior = lda_beta_prior)
#         BlockModel <- with(ret, blocks.pos / (blocks.pos + blocks.neg + 1))
#         MixedMembership <- prop.table(ret$document_expects, 2)
#         colnames(MixedMembership) <- colnames(sociomatrix[[i]])
#         MixedMembership <- MixedMembership[,colnames(MixedMembership) %in% ntid] 
#         temp_res[[i]] <- list(BlockModel = BlockModel,
#                               MixedMembership = MixedMembership)
#         
#         
#       } 
#       else {
#         blk1<- blocks
#         blk2<- ifelse(family2,n.blocks1,n.blocks2)
#         mfd_list <- split(mfd, mfd[,c("(tid)")])
#         mfm_list <- split(mfm, mfm[,c("(tid)")])
#         temp_res[[i]] <- mmsbm(update(formula.dyad, .~1),
#                                formula.monad = ~ 1,
#                                senderID = "(sid)",
#                                receiverID = "(rid)",
#                                nodeID = "(nid)",
#                                timeID = "(tid)",
#                                data.dyad = mfd_list[[i]],
#                                data.monad = mfm_list[[i]][,c("(nid)", "(tid)")],
#                                n.blocks1 = blk1,
#                                n.blocks2 = blk2,
#                                n.hmmstates = 1,
#                                directed = directed,
#                                missing = missing,
#                                mmsbm.control = list(em_iter = em_iter,
#                                                     seed = seed,
#                                                     mu_b = mu_b,
#                                                     var_b = var_b,
#                                                     spectral = spectral,
#                                                     conv_tol = conv_tol,
#                                                     threads = threads,
#                                                     verbose = FALSE))
#         temp_res[[i]]$MixedMembership <- temp_res[[i]]$MixedMembership[,colnames(temp_res[[i]]$MixedMembership) %in% ntid] 
#       }
#     }
#     
#     temp_res <- lapply(split(temp_res, state_init),
#                        function(mods){
#                          target <- t(mods[[1]]$MixedMembership)
#                          rownames(target) <- sapply(strsplit(rownames(target), "@", fixed = TRUE, useBytes = TRUE), function(x)x[1])
#                          res <- lapply(mods,
#                                        function(mod, target_mat = target){
#                                          split_names <- strsplit(colnames(mod$MixedMembership), "@", fixed = TRUE, useBytes = TRUE)
#                                          mod_names <-  sapply(split_names, function(x)x[1])
#                                          mod_time <- split_names[[1]][2]
#                                          shared_nodes <- intersect(mod_names,
#                                                                    rownames(target_mat))
#                                          shared_nodes_mod <- paste(shared_nodes, mod_time, sep="@")
#                                          cost_mat <- mod$MixedMembership[,shared_nodes_mod] %*% target_mat[shared_nodes,]
#                                          perm <- clue::solve_LSAP(t(cost_mat), TRUE)
#                                          mod$MixedMembership <- mod$MixedMembership[perm,]
#                                          mod$BlockModel <- mod$BlockModel[perm, perm]
#                                          return(mod)
#                                        })
#                          return(res)
#                        })
#     block_models <- lapply(temp_res, 
#                            function(mods){
#                              Reduce("+", Map(function(x)x$BlockModel, mods)) / length(mods)
#                            }) 
#     perms_temp <- .findPerm(block_models, use.perms = permute)
#     phis_temp <- Map(function(x)x$MixedMembership, unlist(temp_res, recursive = FALSE))
#     #phi_init_orig <- do.call(cbind,mapply(function(phi,perm){perm %*% phi},
#     #phis_temp, perms_temp[state_init], SIMPLIFY = FALSE))
#     #rownames(phi_init_orig) <- 1:blocks
#     res <- do.call(cbind,mapply(function(phi,perm){perm %*% phi},
#                                 phis_temp, perms_temp[state_init], SIMPLIFY = FALSE))
#     rownames(res) <- 1:blocks
#   }
#   return(res)
# }

#Handle null beta
# .nullBeta <- function(x, state_initial, phi_initial_t, time_id_n){
#   X_state <- split.data.frame(x, state_initial[time_id_n + 1])
#   phi_state <- split.data.frame(t(phi_initial_t), state_initial[time_id_n + 1])
#   res <- mapply(function(X_sub, phi_sub){
#     phi_temp <- .transf(phi_sub)
#     lm.fit(X_sub, log(phi_temp))$coefficients
#   },
#   X_state, phi_state,
#   SIMPLIFY = "array")
#   return(res)
# } 