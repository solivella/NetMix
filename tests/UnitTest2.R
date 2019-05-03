#################################
## Unit test 2:
## One time period, both covariates
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator2.R")
set.seed(831213)
BM <- jitter(matrix(c(5, rep(-5, 3), -5, rep(-5, 3), 1),ncol=3)*-1)
BM_sym <- as.matrix(Matrix::forceSymmetric(BM))
directed <- FALSE
net2 <-  NetSim2(BLK = 3
                ,NODE = 100
                ,STATE = 1
                ,TIME = 1 
                ,DIRECTED = directed
                ,N_PRED=c(1, 1)
                ,B = if(directed) BM else BM_sym
                ,beta_arr = list(array(c(.5, .25,
                                   .5, -1.25, 
                                   0.50, 2.25)*-1,
                                   c(3, 2)))
                ,gamma_vec = c(0.25))

# real_phis2 <- t(do.call(rbind,net2$pi_vecs))
# colnames(real_phis2) <- paste(1:(net2$NODE),1, sep="@")

net2.model <- mmsbm(formula.dyad = Y ~ V1,
                    formula.monad = ~ V2,
                    senderID = "node1",
                    receiverID = "node2",
                    # senderID = "sender",
                    # receiverID = "receiver",
                    nodeID = "node",
                    data.dyad = net2$dyad.data,
                    data.monad = net2$monad.data,
                    n.blocks = net2$BLK,
                    n.hmmstates = 1,
                    directed = net2$DIRECTED,
                    mmsbm.control = list(mu_b = c(2.5,-2.55)*-1,
                                         var_b = c(1, 1),
                                         spectral = TRUE,
                                         verbose = TRUE,
                                         em_iter = 500,
                                         conv_tol = 1e-4
                                         ))
for (i in 1:10) gc()
loss.mat.net2 <- net2.model$MixedMembership %*% net2$pi_vecs#[[1]]
net2_order <- clue::solve_LSAP(t(loss.mat.net2), TRUE)
pred_data_static <- data.frame(Network = rep(c("Network II"), each=net2$NODE*net2$BLK),
                               Group = factor(rep(c(1:3), each = net2$NODE, times = 1)),
                               #Truth = c(c(net2$pi_vecs[[1]])),
                               Truth = c(c(net2$pi_vecs)),
                               Pred = c(c(t(net2.model$MixedMembership[net2_order,]))))
ggplot(pred_data_static, aes(x=Truth, y=Pred, color=Group, pch=Group))+
  scale_color_brewer(palette="Set1") + 
  geom_point(alpha=0.7, size=2.5) + 
  facet_wrap(~Network) + 
  theme_bw() + 
  ylab("Estimate") +
  xlab("True mixed-membership")
##Monad coefs
net2.model$MonadCoef[,net2_order,1]
net2$beta_arr
##Dyad coefs
net2.model$DyadCoef
net2$gamma_mat
## Block model
net2.model$BlockModel[net2_order,net2_order]
net2$B
