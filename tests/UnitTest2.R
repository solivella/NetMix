#################################
## Unit test 2:
## One time period, both covariates
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator.R")
#set.seed(831213)
net2 <-  NetSim(BLK = 3
                ,NODE = 300
                ,STATE = 1
                ,TIME = 1 
                ,DIRECTED = TRUE
                ,N_PRED=1
                ,B_t = matrix(c(5, rep(-5, 3), 5, rep(-5, 3), 5), 
                              ncol=3)*-1
                ,beta_arr = array(c(1.25, 1.25, #BLOCK 1
                                    1.25, -1.25, #BLOCK 2
                                    0.50, 2.25  #BLOCK 3
                                    )*-1,
                                  c(2, 3, 1))
                ,gamma_vec = c(1.5))

real_phis2 <- t(do.call(rbind,net2$pi_vecs))
colnames(real_phis2) <- paste(1:(net2$NODE),1, sep="@")

net2.model <- mmsbm(formula.dyad = Y ~ V1,
                    formula.monad = ~ V1,
                    senderID = "sender",
                    receiverID = "receiver",
                    nodeID = "node",
                    data.dyad = net2$dyad.data,
                    data.monad = net2$monad.data,
                    n.blocks = net2$BLK,
                    n.hmmstates = 1,
                    directed = net2$DIRECTED,
                    mmsbm.control = list(mu_b = c(5,-5)*-1,
                                         var_b = c(1, 1),
                                         spectral = TRUE,
                                         #phi_init_t = real_phis2,
                                         em_iter = 500,
                                         conv_tol = 1e-3
                    ))
loss.mat.net2 <- net2.model$MixedMembership %*% net2$pi_vecs[[1]]
net2_order <- clue::solve_LSAP(t(loss.mat.net2), TRUE)
pred_data_static <- data.frame(Network = rep(c("Network II"), each=net2$NODE*net2$BLK),
                               Group = factor(rep(c(1:3), each = net2$NODE, times = 1)),
                               Truth = c(c(net2$pi_vecs[[1]])),
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
