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
                       function(mat, mat2){
                         res <- matrix(NA, ncol=k, nrow=k)
                         mat <- .reRange(mat)
                         mat2 <- .reRange(mat2)
                         for(i in 1:k){
                           for(j in 1:k){
                             res[i, j] <- sum(.bernKL(mat[i,], mat2[j,]) +
                               sum(.bernKL(mat[,i], mat2[,j])))
                             }
                         }
                         return(res)
                       }, mat2 = B.prime)
    ##Get optimal perm 
    perm.vecs <- lapply(cost.mat,
                        function(mat){
                          clue::solve_LSAP(mat)
                        })

    block.norm.p <- mapply(function(ind, mat){mat[ind,ind]}, perm.vecs, block.norm, SIMPLIFY = FALSE)

    ##Compute risk
    new.risk <- sum(sapply(block.norm.p,
                            function(mat,tar){
                              mat <- .reRange(mat)
                              tar <- .reRange(tar)
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


cluster.mems <- function(fm, t, n=10, demean=FALSE){
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
