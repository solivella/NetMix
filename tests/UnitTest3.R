#################################
## Unit test 3:
## Ten time periods, both covariates
#################################

library(NetMix)
library(tidyverse)
source("tests/NetGenerator.R")
#set.seed(831213)
net3 <-  NetSim(BLK = 3
                ,NODE = 200
                ,STATE = 2
                ,TIME = 10 
                ,DIRECTED = TRUE
                ,N_PRED=2
                ,B_t = matrix(c(5, rep(-5, 3), -5, rep(-5, 3), 1), 
                              ncol=3)*-1
                ,beta_arr = array(c(matrix(c(c(0.25, 0.25, 0.25)*2, #Intercept
                                             0.25, -1.25, 2.25,#Var1
                                           0.25, -1.25, -1.25),#Var2
                                           ncol=3, byrow = TRUE),
                                  matrix(c(c(0.25, 0.25, 0.25)*2, #Intercept
                                           1.25, -1.25, -2.25,#V1
                                         -0.25, -1.25, 2.25),#Var2
                                         ncol=3, byrow = TRUE)),c(3, 3, 2)),
                sVec = rep(c(1,2),c(5,5))
                ,gamma_vec = c(1.5, -1.3))

real_phis3 <- t(do.call(rbind,net3$pi_vecs))
colnames(real_phis3) <- with(subset(net3$monad.data), paste(node,year, sep="@"))

net3.model <- mmsbm(formula.dyad = Y ~ V1 + V2,
      formula.monad = ~ V1 + V2,
      senderID = "sender",
      receiverID = "receiver",
      nodeID = "node",
      timeID = "year",
      data.dyad = net3$dyad.data,
      data.monad = net3$monad.data,
      n.blocks = net3$BLK,
      n.hmmstates = net3$STATE,
      directed = net3$DIRECTED,
      mmsbm.control = list(mu_b = c(5,-5)*-1,
                           var_b = c(1, 1),
                           spectral = TRUE,
                           #phi_init_t = real_phis3,
                           em_iter = 500,
                           conv_tol = 1e-3,
                           threads = 24
                           ))

loss.mat.net3 <- net3.model$MixedMembership %*% do.call(rbind,net3$pi_vecs)
net3_order <- clue::solve_LSAP(t(loss.mat.net3), TRUE)
pred_data_static <- data.frame(Network = rep(c("Network III"), each=net3$NODE*net3$BLK*net3$TIME),
                               Group = factor(rep(c(1:3), each = net3$NODE, times = net3$TIME)),
                               Truth = c(do.call(rbind,net3$pi_vecs)),
                               Pred = c(c(t(net3.model$MixedMembership[net3_order,])))
                               #Year = rep(1:net3$TIME, each=net3$NODE)
                               )
ggplot(pred_data_static, aes(x=Truth, y=Pred, color=Group, pch=Group))+
  scale_color_brewer(palette="Set1") + 
  geom_point(alpha=0.7, size=2.5) + 
  #facet_wrap(~Network) + 
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

##Gof
#gof(net3.model)


