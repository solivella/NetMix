load("tests/MIDdata.Rdata") 
library(NetMix)
set.seed(9999)
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
                                    verbose = TRUE,
                                    em_iter = 3000,
                                    #eta = 1,
                                    threads = 4
                                  )))


