load("tests/MIDdata.Rdata") 
library(NetMix)
#set.seed(9999)
#rand_phis <- prop.table(matrix(runif(length(node_names)*4), c(4, length(node_names))), 2)
#colnames(rand_phis) <- node_names
system.time(fit_onset_pW <- mmsbm(formula.dyad = MID_onset ~ IGOmems_joint + trade_dep_low + 
                                    ally + contiguity + dist +
                                    peaceyrs + spline1 + spline2 + spline3,
                                  formula.monad =  ~ polity + logNMC,
                                  data.dyad = MID_dyad, data.monad = MID_monad,
                                  senderID = "country1", receiverID = "country2", 
                                  nodeID = "country", 
                                  timeID = "year",
                                  n.groups = 4, 
                                  n.hmmstates = 3,
                                  directed=FALSE,
                                  mmsbm.control = list(
                                    init = "spectral",
                                    verbose = TRUE,
                                    em_iter = 3000,
                                    conv_tol = 1e-6,
                                    #phi_init_t = rand_phis,
                                    #eta = 1,
                                    threads = 4
                                  )))


