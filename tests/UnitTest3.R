#################################
## Unit test 3:
## Ten time periods, both covariates
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator2.R")
rseed <- as.numeric(paste(sample(1:5), collapse=""))
set.seed(rseed)
BM <- jitter(matrix(c(5, rep(-5, 3), -5, rep(-5, 3), 1),ncol=3)*-1)
BM_sym <- as.matrix(Matrix::forceSymmetric(BM))
directed <- FALSE
net3 <-  NetSim2(BLK = 3
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

real_phis3 <- t(net3$pi_vecs)
colnames(real_phis3) <- with(subset(net3$monad.data), paste(node,time, sep="@"))

# net3$dyad.data2 <- net3$dyad.data[net3$dyad.data$node1 != net3$dyad.data$node2,]
# net3$dyad.data3 <- net3$dyad.data
# net3$dyad.data3$Y[net3$dyad.data3$node1 != net3$dyad.data3$node2] <- 0


net3.model <- mmsbm(formula.dyad = Y ~  V1,
                    formula.monad = ~ V2, ## Change this to V2 of using new simulaton code
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
                                         seed=123,
                                         alpha = 1.0,
                                         #b_init_t = t(net3$B),
                                         #beta_init = sapply(net3$beta_arr, function(x)t(x), simplify = "array"),
                                         #gamma_init = net3$gamma_mat,
                                         #phi_init_t = real_phis3,
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
#pdf("Net3Test.pdf")
ggplot(pred_data_static, aes(x=Truth, y=Pred, color=Group, pch=Group))+
  scale_color_brewer(palette="Set1") + 
  geom_point(size=2.5, alpha=0.7) + 
  facet_wrap(~Year) + 
  theme_bw() + 
  ylab("Estimate") +
  xlab("True mixed-membership")
#dev.off()
##Monad coefs
# net3.model$MonadCoef[,net3_order,]
# net3$beta_arr
# ##Dyad coefs
# net3.model$DyadCoef
# net3$gamma_mat
# ## Block model
net3.model$BlockModel[net3_order,net3_order]
net3$B
# ## Kappa
round(t(net3.model$Kappa),2)
net3$sVec
# 
# 
