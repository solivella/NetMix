load("tests/MIDdata.Rdata") 
library(NetMix)
set.seed(9999)
fit_onset_pW <- mmsbm(formula.dyad = MID_onset ~ IGOmems_joint + trade_dep_low + 
                        ally + contiguity + dist +
                        peaceyrs + spline1 + spline2 + spline3,
                      formula.monad =  ~ polity + logNMC + gdp_growth_67_a + region,
                      data.dyad = MID_dyad, data.monad = MID_monad,
                      senderID = "country1", receiverID = "country2", 
                      nodeID = "country", 
                      timeID = "year",
                      n.groups = 4, 
                      n.hmmstates = 2,
                      directed=FALSE,
                      mmsbm.control = list(var_b = c(9,9),
                                           var_xi = 10,
                                           var_gamma = 9,
                                           mu_b = c(-5, 5),
                                           verbose = TRUE,
                                           em_iter = 3000,
                                           eta = .1,
                                           threads = parallel::detectCores()
                      ))
