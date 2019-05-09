#################################
## Unit test 4:
## Ten time periods, both covariates, missing nodes
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator2.R")
rseed <- as.numeric(paste(sample(1:5), collapse=""))
set.seed(rseed)
BM <- jitter(matrix(c(5, rep(-5, 3), -5, rep(-5, 3), 1),ncol=3)*-1)
BM_sym <- as.matrix(Matrix::forceSymmetric(BM))
directed <- FALSE
net4 <-  NetSim2(BLK = 3
                 ,NODE = 100
                 ,STATE = 2
                 ,TIME = 10
                 ,DIRECTED = directed
                 ,N_PRED=c(1, 1)
                 ,B =  if(directed) BM else BM_sym
                 ,beta_arr = list(array(c(-1.25, -.75,
                                          -1.25, 0.375,
                                          -1.25, -4.625),
                                        c(3, 2)),
                                  array(c(1.25, 1.375,
                                          -1.25, 0.625,
                                          1.25, 4.625),
                                        c(3, 2)))
                 ,sVec = rep(c(1,2),c(5,5))
                 ,gamma_vec = c(0.05)
)


real_phis4 <- t(net4$pi_vecs)
colnames(real_phis4) <- with(subset(net4$monad.data), paste(node,time, sep="@"))

node.rem1 <- sample(unique(net4$monad.data$node), 15)
node.rem2 <- sample(unique(net4$monad.data$node), 10)

#net4$monad.data2 <- net4$monad.data
#net4$dyad.data2 <- net4$dyad.data
net4$monad.data2 <- net4$monad.data[-c(which(net4$monad.data$node %in% node.rem1 & net4$monad.data$time <= 3),
                                       which(net4$monad.data$node %in% node.rem2 & net4$monad.data$time <= 7)),]
net4$dyad.data2 <- net4$dyad.data[-c(which((net4$dyad.data$node1 %in% node.rem1 | net4$dyad.data$node2 %in% node.rem1) & 
                                       net4$dyad.data$time <= 3),
                                     which((net4$dyad.data$node1 %in% node.rem2 | net4$dyad.data$node2 %in% node.rem2) & 
                                             net4$dyad.data$time <= 7)),]

net4.model <- mmsbm(formula.dyad = Y ~  V1,
                    formula.monad = ~ V2, 
                    senderID = "node1",
                    receiverID = "node2",
                    nodeID = "node",
                    timeID = "time",
                    data.dyad = net4$dyad.data2,
                    data.monad = net4$monad.data2,
                    n.blocks = net4$BLK,
                    n.hmmstates = net4$STATE,
                    directed = net4$DIRECTED,
                    mmsbm.control = list(mu_b = c(5,-5)*-1,
                                         var_b = c(1, 1),
                                         spectral = TRUE,
                                         seed=123,
                                         alpha = 1.0,
                                         #b_init_t = t(net4$B),
                                         #beta_init = sapply(net4$beta_arr, function(x)t(x), simplify = "array"),
                                         #gamma_init = net4$gamma_mat,
                                         #phi_init_t = real_phis3,
                                         #kappa_init_t = t(model.matrix(~-1+as.factor(net4$sVec))),
                                         verbose=TRUE,
                                         em_iter = 300,
                                         conv_tol = 1e-4,
                                         threads = 4
                    ))

#net4$pi_vecs2 <- net4$pi_vecs
net4$pi_vecs2 <- net4$pi_vecs[paste(net4$monad.data2$node, net4$monad.data2$time, sep="_"),]

loss.mat.net4 <- net4.model$MixedMembership %*% net4$pi_vecs2

net4_order <-  clue::solve_LSAP(t(loss.mat.net4), TRUE)

pred_data_static <- data.frame(Network = rep(c("Network III"), each=nrow(net4$monad.data2)*net4$BLK),
                               Group = factor(rep(c(1:net4$BLK), each = nrow(net4$monad.data2))),
                               Year = rep(rep(1:10, times=table(net4$monad.data2$time)), times=net4$BLK),
                               Truth = c(net4$pi_vecs2),
                               Pred = c(t(net4.model$MixedMembership[net4_order,])))
#pdf("net4Test.pdf")
ggplot(pred_data_static, aes(x=Truth, y=Pred, color=Group, pch=Group))+
  scale_color_brewer(palette="Set1") + 
  geom_point(size=2.5, alpha=0.7) + 
  facet_wrap(~Year) + 
  theme_bw() + 
  ylab("Estimate") +
  xlab("True mixed-membership")
#dev.off()
##Monad coefs
# net4.model$MonadCoef[,net4_order,]
# net4$beta_arr
# ##Dyad coefs
# net4.model$DyadCoef
# net4$gamma_mat
# ## Block model
net4.model$BlockModel[net4_order,net4_order]
net4$B
# ## Kappa
round(t(net4.model$Kappa),2)
net4$sVec
# 
# 
