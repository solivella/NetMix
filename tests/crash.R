library(NetMix)
load("tests/fit_onset.Rdata")
fit_onset <- MID_onset_6_2o


btest <- boot.mmsbm(fit_onset, 5)
