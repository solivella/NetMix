#################################
## Unit test 3:
## Ten time periods, both covariates
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator2.R")
source("tests/NetGenerator.R")
set.seed(831213)
net3 <-  NetSim2(BLK = 3
                 ,NODE = 100
                 ,STATE = 2
                 ,TIME = 10
                 ,DIRECTED = FALSE
                 ,N_PRED=c(1, 1)
                 ,B = jitter(matrix(c(-5, rep(5, 3), -3, rep(5, 3), -1),
                                    ncol=3))
                 ,beta_arr = list(array(c(.5, .25,
                                          .5, -1.25,
                                          0.50, 2.25)*-1,
                                        c(3, 2)),
                                  array(c(1.5, .55,
                                          .5, -1.25,
                                          0.50, -1.25)*-1,
                                        c(3, 2))),
                 sVec = rep(c(1,2),c(5,5))
                 ,gamma_vec = c(1.5))

net3 <-  NetSim(BLK = 3
                ,NODE = 100
                ,STATE = 2
                ,TIME = 10
                ,DIRECTED = TRUE
                ,N_PRED=c(1, 1)
                ,B = jitter(matrix(c(-5, rep(5, 3), -3, rep(5, 3), -1),
                                   ncol=3))
                ,beta_arr = list(array(c(.5, .25,
                                         .5, -1.25,
                                         0.50, 2.25)*-1,
                                       c(3, 2)),
                                 array(c(1.5, .55,
                                         .5, -1.25,
                                         0.50, -1.25)*-1,
                                       c(3, 2))),
                sVec = rep(c(1,2),c(5,5))
                ,gamma_vec = c(1.5))

real_phis3 <- t(net3$pi_vecs)
colnames(real_phis3) <- with(subset(net3$monad.data), paste(node,time, sep="@"))

net3.model <- mmsbm(formula.dyad = Y ~  V1,
                    formula.monad = ~ V1, ## Change this to V2 of using new simulaton code
                    senderID = "node1",
                    receiverID = "node2",
                    #senderID = "sender",
                    #receiverID = "receiver",
                    nodeID = "node",
                    timeID = "time",
                    data.dyad = net3$dyad.data,
                    data.monad = net3$monad.data,
                    n.blocks = net3$BLK,
                    n.hmmstates = net3$STATE,
                    directed = net3$DIRECTED,
                    mmsbm.control = list(mu_b = c(5,-5)*-1,
                                         var_b = c(1, 1),
                                         spectral = TRUE,
                                         b_init_t = net3$B,
                                         gamma_init = net3$gamma_mat,
                                         phi_init_t = real_phis3,#[,(x*net3$NODE + 1):(x*net3$NODE + net3$NODE)],
                                         #kappa_init_t = t(model.matrix(~-1+as.factor(net3$sVec))),
                                         verbose=TRUE,
                                         em_iter = 300,
                                         conv_tol = 1e-4,
                                         threads = 4
                    ))

loss.mat.net3 <- net3.model$MixedMembership %*% net3$pi_vecs

net3_order <-  clue::solve_LSAP(t(loss.mat.net3), TRUE)




pred_data_static <- data.frame(Network = rep(c("Network III"), each=net3$NODE*net3$BLK*net3$TIME),
                               Group = factor(rep(c(1:net3$BLK), each = net3$NODE, times = net3$TIME)),
                               Year = rep(1:10, each=net3$NODE, times=net3$BLK),
                               Truth = c(net3$pi_vecs),
                               Pred = c(t(net3.model$MixedMembership[net3_order,])))

ggplot(pred_data_static, aes(x=Truth, y=Pred, color=Group, pch=Group,alpha=Year))+
  scale_color_brewer(palette="Set1") + 
  geom_point(size=2.5) + 
  facet_wrap(~Network) + 
  theme_bw() + 
  ylab("Estimate") +
  xlab("True mixed-membership")
##Monad coefs
net3.model$MonadCoef[,net3_order,]
net3$beta_arr
##Dyad coefs
net3.model$DyadCoef
net3$gamma_mat
## Block model
net3.model$BlockModel[net3_order,net3_order]
net3$B
## Kappa
t(net3.model$Kappa)
net3$sVec


