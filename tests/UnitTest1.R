#################################
## Unit test 1:
## One time period, no covariates
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator2.R")
set.seed(831213)
BM <- jitter(matrix(c(5, rep(-5, 3), -5, rep(-5, 3), 1),ncol=3)*-1)
BM_sym <- as.matrix(Matrix::forceSymmetric(BM))
directed <- TRUE
net1 <- NetSim2(BLK = 3
                 ,NODE = 50
                 ,STATE = 1
                 ,TIME = 1
                 ,N_PRED = c(0, 1)
                ,DIRECTED = directed
                ,B = if(directed) BM else BM_sym
                ,beta_arr = list(array(c(1.25, 1.25, 1.25)*-1, c(3,1)))#array(c(1.25, 1.25, 1.25)*-1,c(1, 3, 1))
                ,gamma_vec = c(-0.25)
)

#real_phis1 <- t(do.call(rbind,net1$pi_vecs))
#colnames(real_phis1) <- paste(1:(net1$NODE),1, sep="@")

net1.model <- mmsbm(formula.dyad = Y ~ V1,
                    senderID = "node1",
                    receiverID = "node2",
                    #senderID = "sender",
                    #receiverID = "receiver",
                    #timeID = "time",
                    data.dyad = net1$dyad.data,
                    n.blocks = net1$BLK,
                    n.hmmstates = net1$STATE,
                    directed = net1$DIRECTED,
                    mmsbm.control = list(mu_b = c(5,-5)*-1,
                                         var_b = c(1, 1),
                                         spectral = TRUE,
                                         verbose=TRUE,
                                         em_iter = 5000,
                                         conv_tol = 1e-4
                    ))
#loss.mat.net1 <- net1.model$MixedMembership %*% net1$pi_vecs[[1]]
loss.mat.net1 <- net1.model$MixedMembership %*% net1$pi_vecs
net1_order <- clue::solve_LSAP(t(loss.mat.net1), TRUE)
pred_data_static <- data.frame(Network = rep(c("Network I"), each=net1$NODE*net1$BLK),
                               Group = factor(rep(c(1:3), each = net1$NODE, times = 1)),
                              #Truth = c(c(net1$pi_vecs[[1]])),
                              Truth = c(c(net1$pi_vecs)),
                               Pred = c(c(t(net1.model$MixedMembership[net1_order,]))))
ggplot(pred_data_static, aes(x=Truth, y=Pred, color=Group, pch=Group))+
  scale_color_brewer(palette="Set1") + 
  geom_point(alpha=0.7, size=2.5) + 
  facet_wrap(~Network) + 
  theme_bw() + 
  ylab("Estimate") +
  xlab("True mixed-membership")
##Monad coefs
net1.model$MonadCoef[,net1_order,1]
net1$beta_arr
##Dyad coefs
net1.model$DyadCoef
net1$gamma_mat
## Block model
net1.model$BlockModel[net1_order,net1_order]
net1$B


