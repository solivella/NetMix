library(NetMix)
source("Extra/NetGenerator.R")
set.seed(9999)
sim <- NetSim(BLK=3, NODE=30, STATE=2, TIME=10, DIRECTED=TRUE, N_PRED=3,
              beta_arr = array(rnorm(12), dim=c(4,3,2)),
              gamma_vec = rnorm(3, 2, 0.1))

set.seed(9999)
f <-  mmsbm(formula.dyad = Y ~ 1,
            formula.monad = ~ V1 + V2 + V3,
            data.dyad = sim$dyad.data, 
            data.monad = sim$monad.data,
            senderID="sender", receiverID="receiver",
            nodeID="node", timeID="year",
            n.groups=sim$BLK, n.hmmstates=sim$STATE,
            directed=TRUE,
            mmsbm.control = list(verbose = TRUE,
                                 threads = parallel::detectCores()))

