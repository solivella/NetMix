#########################
## Misc. Helper functions
#########################

.modSum <- function(x, y){
  replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
}

.bernKL <- function(x,y){
  x * log(x/y) + (1-x)*log((1-x)/(1-y))
}

.findPerm <- function(block.list, max.iter, target.mat = NULL){
  n <- length(block.list)
  k <- ncol(block.list[[1]])
  block.norm <- block.list
  block.norm.p <- block.norm
  if(!is.null(target.mat)){
    target.mat <- prop.table(target.mat, 1)
  }
  iter <- 0
  old.risk <- 1e6
  chg <- 1
  while((chg > 1e-6) && (iter < max.iter)){
    ## Step 1
    if(is.null(target.mat)){
      B.prime <- Reduce("+",block.norm.p)/n
    } else {
      B.prime <- target.mat
    }
    ## Step 2
    ## Compute loss matrix (KL)
    cost.mat <- lapply(block.norm,
                       function(mat){
                         res <- matrix(NA, ncol=k, nrow=k)
                         if(any(mat<1e-12)){mat[mat < 1e-12] <- 1e-12}
                         for(i in 1:k){
                           for(j in 1:k){
                             res[i, j] <- sum(.bernKL(mat[i,], B.prime[j,]) +
                               sum(.bernKL(mat[,i], B.prime[,j])))
                             }
                         }
                         return(res)
                       })
    ##Get optimal perm 
    perm.vecs <- lapply(cost.mat,
                        function(mat){
                          clue::solve_LSAP(mat)
                        })

    block.norm.p <- mapply(function(ind, mat){mat[ind,ind]}, perm.vecs, block.norm, SIMPLIFY = FALSE)

    ##Compute risk
    new.risk <- sum(sapply(block.norm.p,
                            function(mat,tar){
                              if(any(mat<1e-12)){mat[mat < 1e-12] <- 1e-12}
                              sum(.bernKL(mat, tar))
                            },
                            tar = B.prime))

    chg <- abs(new.risk - old.risk)
    old.risk <- new.risk
    iter <- iter + 1
  }
  perm.vecs
}

.transf <- function(mat){
  (mat * (nrow(mat) - 1) + 1/ncol(mat))/nrow(mat)
}
.mapID <- function(uvec, vec){
  temp_map <- seq_along(uvec)
  names(temp_map) <- as.character(uvec)
  temp_map[as.character(vec)]
}

## Define functions for displaying model results

summary.mmsbm <- function(fm){
  monad <- fm$MonadCoef
  states <- rowMeans(fm$Kappa)
  summ <- list(nrow(fm$dyadic.data), ncol(fm$BlockModel), 
               rowMeans(fm$MixedMembership),
               exp(fm$BlockModel) / (1 + exp(fm$BlockModel)), 
               fm$DyadCoef, monad, states)
  names(summ) <- c("N", "Number of Clusters", "Percent of Observations in Each Cluster",
                   "Edge Formation Probabilities", "Dyadic Coefficients", "Monadic Coefficients",
                   "Markov State Probabilities")
  print(summ)
}



plot.mmsbm <- function(fm, directed=FALSE){ # network graph showing B-matrix
  mode <- ifelse(directed, "directed", "undirected")
  require("igraph", quietly=TRUE)
  block.G <- graph.adjacency(exp(fm$BlockModel) / (1 + exp(fm$BlockModel)), mode=mode, weighted=TRUE)
  linMap <- function(x, from, to){(x - min(x)) / max(x - min(x)) * (to - from) + from}
  e.weight <- E(block.G)$weight*100
  if(any(e.weight < 0.1)){
    e.weight <- e.weight*200
  }
  if(any(e.weight > 20)){
    e.weight[e.weight > 20] <- 20}
  v.size <- rowMeans(fm$MixedMembership)*150
  plot(block.G, main = "Edge Formation across Clusters",
       edge.width=e.weight, vertex.size=v.size,
       layout = layout_in_circle)
}



cluster.mems <- function(fm, t, n=10, demean=FALSE){
  Mem <- fm$MixedMembership[,fm$monadic.data[,"(tid)"] %in% t]
  Nodes <- fm$monadic.data[,"(nid)"][fm$monadic.data[,"(tid)"] %in% t]
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


cluster.time.node2 <- function(fm){
  for(i in 1:nrow(fm$MixedMembership)){
    plot(unique(fm$monadic.data[,"(tid)"]), 1:nrow(fm$Kappa), cex=0,
         ylim=c(0,1), xlim=c(min(fm$monadic.data[,"(tid)"]), max(fm$monadic.data[,"(tid)"])+1),
         xlab="Time", ylab="Group Membership", main=paste("Membership in Group", i))
    for(j in unique(fm$monadic.data[,"(nid)"])){
      lines(unique(fm$monadic.data[,"(tid)"]), sapply(unique(fm$monadic.data[,"(tid)"]), function(x){
        fm$MixedMembership[i,fm$monadic.data[,"(tid)"]==x & fm$monadic.data[,"(nid)"]==j]}), col=j)
      text(as.character(j), x=max(unique(fm$monadic.data[,"(tid)"]))+1, col=j,
           y=fm$MixedMembership[i,fm$monadic.data[,"(tid)"]==max(unique(fm$monadic.data[,"(tid)"])) & fm$monadic.data[,"(nid)"]==j])
      
    }
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



radar.plot <- function(fm){
  if(nrow(fm$MixedMembership)>2){
    require(fmsb, quietly=TRUE)
    for(i in unique(fm$monadic.data[,"(nid)"])){
      radarchart(rbind(rep(1, 3), rep(0, 3), as.data.frame(t(fm$MixedMembership[,fm$monadic.data[,"(nid)"]==i]))),
                 plty=1, pcol=rgb(0,0,0, alpha=0.1), pfcol=rgb(0,0,0, alpha=0.1),
                 vlabels=paste("Group", 1:nrow(fm$MixedMembership)), title=as.character(i))
    }
  }
  if(nrow(fm$MixedMembership)<=2){print("Number of Groups must be 3 or more")}
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


cluster.top <- function(fm, num){
  top <- sapply(1:fm$n_blocks, function(x){
    ref <- order(fm$MixedMembership[x,], decreasing=TRUE)[1:num]
    paste(fm$monadic.data$`(nid)`[ref], fm$monadic.data$`(tid)`[ref])
  })
  print(paste("Top", num, "node-years by group membership"))
  print(top)
}


est.edgeplot <- function(fm){
  pi1 <- fm$MixedMembership[,match(fm$dyadic.data[,"(sid)"],fm$monadic.data[,"(nid)"])]
  pi2 <- fm$MixedMembership[,match(fm$dyadic.data[,"(rid)"],fm$monadic.data[,"(nid)"])]
  est.ties <- rep(0, ncol(pi1))
  for(a in 1:nrow(pi1)){ # vectorize this
    for(b in 1:nrow(pi2)){
      est.ties <- est.ties + (pi1[a,]*pi2[b,]*fm$BlockModel[a,b]) # NAs in pi?
    }
  }
  est.ties <- est.ties  + t(fm$DyadCoef)%*%t(fm$dyadic.data[,5:ncol(fm$dyadic.data)])
  estimated <- exp(est.ties) / (1 + exp(est.ties))
  obs.node <- sapply(unique(fm$dyadic.data[,"(sid)"]), function(x){
    sum(fm$dyadic.data[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x,4])})
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
  pi1 <- fm$MixedMembership[,match(fm$dyadic.data[,"(sid)"],fm$monadic.data[,"(nid)"])]
  pi2 <- fm$MixedMembership[,match(fm$dyadic.data[,"(rid)"],fm$monadic.data[,"(nid)"])]
  est.ties <- rep(0, ncol(pi1))
  for(a in 1:nrow(pi1)){
    for(b in 1:nrow(pi2)){
      est.ties <- est.ties + (pi1[a,]*pi2[b,]*fm$BlockModel[a,b]) + t(fm$DyadCoef)%*%t(fm$dyadic.data[,5:ncol(fm$dyadic.data)])
    }
  }
  estimated <- exp(est.ties) / (1 + exp(est.ties))
  obs.node <- sapply(unique(fm$dyadic.data[,"(sid)"]), function(x){
    sum(fm$dyadic.data[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x,4])})
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



est.pi <- function(fm, monad=fm$monadic.data){
  pi.states <- lapply(1:nrow(fm$Kappa), function(m){
    ifelse(nrow(fm$MonadCoef) > 1, 
           return((t(fm$MonadCoef[,,m]) %*% t(cbind(rep(1,nrow(monad)), monad[,1:(ncol(monad)-2)]))) *
                    fm$Kappa[m, as.numeric(factor(monad[,"(tid)"]))]),
           return((as.matrix(fm$MonadCoef[,,m]) %*% t(rep(1,nrow(monad)))) * 
                    fm$Kappa[m, as.numeric(factor(monad[,"(tid)"]))])
    )
  })
  return(Reduce("+", pi.states))
}

predict.mmsbm <- function(fm, dyad=fm$dyadic.data, monad=fm$monadic.data, 
                          out.sample=FALSE, posterior.pi=FALSE,
                          type="probability"){ 
  if(!out.sample){
    ifelse(posterior.pi, pis <- fm$MixedMembership, pis <- est.pi(fm, monad))
    sid <- ifelse(identical(fm$dyadic.data, dyad), "(sid)", fm$call$senderID)
    rid <- ifelse(identical(fm$dyadic.data, dyad), "(rid)", fm$call$receiverID)
    nid <- ifelse(identical(fm$dyadic.data, dyad), "(nid)", fm$call$nodeID)
    pi1 <- pis[,match(dyad[,sid],monad[,nid])]
    pi2 <- pis[,match(dyad[,rid],monad[,nid])]
  }
  
  if(out.sample){
    tid <- ifelse(all(colnames(fm$monadic.data) %in% colnames(monad)), "(tid)", fm$call$timeID)
    kappa.last <- fm$Kappa[,ncol(fm$Kappa)]
    kappa.new <- list()
    kappa.new[[1]] <- fm$TransitionKernel %*% kappa.last
    if(length(unique(dyad$year))>1){
      for(i in 2:length(unique(dyad$year))){
        kappa.new[[i]] <- fm$TransitionKernel %*% kappa.new[[i-1]]
      }
    }
    #kappa.new[[1]] <- kappa.new[[2]] <- c(1, 0, 0)
    new.pis <- list()
    for(i in 1:fm$n_states){
      ifelse(!is.null(fm$call$formula.monad),
             monad.vars <- cbind(rep(1, nrow(monad)), get_all_vars(fm$call$formula.monad, monad)),
             monad.vars <- rep(1, nrow(fm$monadic.data)))
      new.pis[[i]] <- as.matrix(monad.vars) %*% fm$MonadCoef[,,i]
    }
    mean.pis <- new.pis[[1]]
    for(j in 1:length(unique(monad[,tid]))){
      kappa.rel <- kappa.new[[j]]
      mean.pi.sub <- mean.pis[monad[,tid]==unique(monad[,tid])[j],]
      for(q in 1:nrow(mean.pi.sub)){
        pi.mean <- rep(0, ncol(mean.pi.sub))
        for(v in 1:fm$n_states){
          pi.mean <- pi.mean + new.pis[[v]][q,] * kappa.new[[j]][v]
        }
        mean.pi.sub[q,] <- pi.mean
      }
      mean.pis[monad[,tid]==unique(monad[,tid])[j],] <- mean.pi.sub
    }
    pi1 <- t(mean.pis[match(dyad[,fm$call$senderID], monad[,fm$call$nodeID]),])
    pi2 <- t(mean.pis[match(dyad[,fm$call$receiverID], monad[,fm$call$nodeID]),])
  }
  
  est.ties <- rep(0, ncol(pi1))
  for(a in 1:nrow(pi1)){ 
    for(b in 1:nrow(pi2)){
      est.ties <- est.ties + (pi1[a,]*pi2[b,]*fm$BlockModel[a,b])
    }
  }
  est.ties <- est.ties  + t(fm$DyadCoef)%*%t(get_all_vars(fm$call$formula.dyad, dyad)[,-1])
  if(type=="probability"){
    return(exp(est.ties) / (1 + exp(est.ties)))}
  if(type=="response"){
    return(rbinom(length(est.ties), 1, prob=exp(est.ties) / (1 + exp(est.ties))))
  }
}


covFX <- function(fm, cov, shift){
  predict.ties <- predict(fm)
  monadic.data2 <- fm$monadic.data
  monadic.data2[,cov] <- fm$monadic.data[,cov] + 1
  predict.ties2 <- predict(fm, monad=monadic.data2)
  FX <- list(mean(predict.ties2 - predict.ties), #avg
             tapply(predict.ties2-predict.ties, fm$dyadic.data[,3], mean), #time
             sapply(unique(fm$monadic.data[,"(nid)"]), function(x){ #node
               mean((predict.ties2-predict.ties)[fm$dyadic.data[,"(sid)"]==x | fm$dyadic.data[,"(rid)"]==x])}),
             tapply(predict.ties2-predict.ties, paste(fm$dyadic.data[,"(sid)"], fm$dyadic.data[,"(rid)"], sep="_"), mean),#dyad
             predict.ties2 - predict.ties) #dyad-time
  names(FX[[5]]) <- paste(fm$dyadic.data[,"(sid)"], fm$dyadic.data[,"(rid)"], sep="_")
  names(FX) <- c(paste("Overall Avg. Effect of", cov), paste("Avg. Effect of", cov, "by Time"),
                 paste("Avg. Effect of", cov, "by Node"), paste("Avg. Effect of", cov, "by Dyad"),
                 paste("Effect of", cov, "by Dyad-Time"))
  return(FX)
}


plot.FX <- function(FX, fm){
  cov <- strsplit(names(FX)[1], " ")[[1]][5]
  ymax <- max(hist(FX[[5]])[["counts"]])
  hist(FX[[5]], main=paste("Distribution of Marginal Effects:", strsplit(names(FX)[1], " ")[[1]][5]),
       xlab=paste("Effect of", cov, "on Pr(Edge Formation)"))
  lines(x=c(FX[[1]], FX[[1]]), y=c(0,ymax*1.05), col="red", lwd=2)
  text(x=FX[[1]], y=ymax, paste("Avg. Effect =", round(FX[[1]],4)), col="red", pos=4)
  
  plot(unique(fm$dyadic.data[,"(tid)"]), tapply(FX[[5]], fm$dyadic.data[,"(tid)"], mean), type="o",
       xlab="Time", ylab=paste("Effect of", cov, "on Pr(Edge Formation)"), main="Marginal Effect over Time")
  
  dyads <- names(FX[[4]])
  totals <- cbind(as.character(unlist(lapply(strsplit(dyads, "_"), '[[', 1))),
                  as.character(unlist(lapply(strsplit(dyads, "_"), '[[', 2))),
                  as.vector(tapply(FX[[5]], names(FX[[5]]), mean)))
}





boot.mmsbm <- function(fm, Y, cov){  ## add multiple core function here (for each loop)
  results <- list()
  results.cov <- list()
  if(ncol(cov) != ncol(fm$coef)-1){cov <- cov[,match(colnames(fm$coef)[2:ncol(fm$coef)], colnames(cov))]}
  orig.order <- colSums(fm$pi) # max, min, middle (will need to adjust this)
  for(i in 1:1000){
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
