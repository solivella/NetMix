library(tidyverse); library(NetMix); library(countrycode)

## Prepare the data

#SIGO <- qs::qread("~/Dropbox/GOV/Research/Network/SIGO3.qs")
SIGO_clean <- SIGO %>% 
  filter(year >= 1950 & year != 1964) %>% 
  mutate(year = ifelse(year < 1964, year - 1949, year - 1950)) %>% 
  mutate(year = as.character(year), ccode = as.character(ccode)) %>% 
  mutate(East_Asia_and_Pacific = ifelse(state_region == "East Asia and Pacific", 1, 0),
         Europe_and_Central_Asia = ifelse(state_region == "Europe and Central Asia", 1, 0),
         Latin_America_and_the_Caribbean = ifelse(state_region == "Latin America and the Caribbean", 1, 0),
         Middle_East_and_North_Africa = ifelse(state_region == "Middle East and North Africa", 1, 0),
         North_America = ifelse(state_region == "North America", 1, 0),
         South_Asia = ifelse(state_region == "South Asia", 1, 0),
         Sub_Saharan_Africa = ifelse(state_region == "Sub-Saharan Africa", 1, 0))


state_monad <- SIGO_clean %>% 
  dplyr::select(c("ccode", "year", "UN_IP", "Polity", "GDPpc",
                  "East_Asia_and_Pacific", "Europe_and_Central_Asia",
                  "Latin_America_and_the_Caribbean", "Middle_East_and_North_Africa",
                  "North_America", "South_Asia", "Sub_Saharan_Africa", "trade.openness2")) %>% 
  distinct(ccode, year, .keep_all = T) %>% 
  rename(VarS1 = UN_IP, VarS2 = Polity, VarS3 = GDPpc, VarS4 = East_Asia_and_Pacific,
         VarS5 = Europe_and_Central_Asia, VarS6 = trade.openness2, 
         id = ccode)

IGO_monad <- SIGO_clean %>% 
  dplyr::select(c("IGO", "year", "Security_IGO", "leadstate.dynamic_authoritarian", "regional_org", 
                  "Econ_IGO", "Environ_IGO", "salient", "IGO_mems_lag_form")) %>% 
  distinct(IGO, year, .keep_all = T) %>% 
  rename(VarB1 = regional_org, VarB2 = leadstate.dynamic_authoritarian,
         VarB3 = Security_IGO, VarB4 = Econ_IGO, VarB5 = Environ_IGO, VarB6 = IGO_mems_lag_form,
         id = IGO)

SIGO_dyad <- SIGO_clean %>% 
  dplyr::select(c("ccode","IGO", "year","Member", "alliances_avgmembers")) %>% 
  distinct(ccode, IGO, year, .keep_all = T) %>% 
  rename(var1 = alliances_avgmembers, Y = Member, id1 = ccode, id2 = IGO)

netSim <- list(df_dyad_1=SIGO_dyad,
               df_monad_B=IGO_monad,df_monad_S=state_monad)

## Fitting the model
beta_mu_array_s <- array(c(0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0),
                         c(7, 3, 2))

beta_mu_array_b <- array(c(0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0),
                         c(7, 3, 2))

beta_var_array_s <- array(c(0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1), 
                          c(7, 3, 2)) #ncov, group, state

beta_var_array_b <- array(c(0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1,
                            0.0001, 1, 1, 1, 1, 1, 1), 
                          c(7, 3, 2)) #ncov, group, state

SIGO_dynbi <- mmsbm(formula.dyad = Y~1, #var1,
                    formula.monad = list(~VarS1 + VarS2 + VarS3 + VarS4 + VarS5 + VarS6, 
                                         ~VarB1 + VarB2 + VarB3 + VarB4 + VarB5 + VarB6),
                   timeID = "year",
                   senderID = "id1",
                   receiverID = "id2",
                   nodeID = list("id","id"),
                   bipartite= TRUE,
                   data.dyad = netSim[["df_dyad_1"]],
                   data.monad = list(netSim[["df_monad_S"]],netSim[["df_monad_B"]]),
                   n.blocks = c(3,3), n.hmmstates = 2, 
                   #moretimes=TRUE, realign=TRUE,fp5times = TRUE,
                   mmsbm.control = list(verbose = TRUE,
                                        threads=1,
                                        svi = FALSE,
                                        vi_iter = 8000,
                                        batch_size = 1.0,
                                        conv_tol = 1e-4,
                                        mu_beta = list(beta_mu_array_s,
                                                        beta_mu_array_b),
                                        var_beta = list(beta_var_array_s,
                                                        beta_var_array_b),
                                        hessian = FALSE))
summary.mmsbmB(SIGO_dynbi)
plot.mmsbmB(SIGO_dynbi, type = "group")  
plot.mmsbmB(SIGO_dynbi, type = "membership_1")  
plot.mmsbmB(SIGO_dynbi, type = "membership_2")  
plot.mmsbmB(SIGO_dynbi, type = "hmm") 





## Look at individual nodes in both families

state_mem <- SIGO_dynbi$MixedMembership1 %>%
  as.data.frame() %>% 
  rownames_to_column(var = "group") %>%
  pivot_longer(cols = -group, names_to = "key", values_to = "value") %>%
  separate(key, into = c("state", "year"), sep = "@") %>% 
  mutate(year = as.numeric(year), state = as.numeric(state),
         year = as.numeric(year) + 1949)

state_top <- state_mem %>% 
  group_by(group, state) %>% 
  summarise_at("value", mean, na.rm = TRUE) %>% 
  mutate(state = countrycode(state, "cown", "country.name"),
         group_rank = rank(-value)) %>% 
  filter(group_rank <= 50)

io_mem <- SIGO_dynbi$MixedMembership2 %>%
  as.data.frame() %>% 
  rownames_to_column(var = "group") %>%
  pivot_longer(cols = -group, names_to = "key", values_to = "value") %>%
  separate(key, into = c("IO", "year"), sep = "@") %>% 
  mutate(year = as.numeric(year) + 1949)
  
io_top <- io_mem %>% 
  group_by(group, IO) %>% 
  summarise_at("value", mean, na.rm = TRUE) %>% 
  mutate(group_rank = rank(-value)) %>% 
  filter(group_rank <= 20)

IO_names <- SIGO %>% select(c(IGO, IGOname)) %>% distinct() %>% rename(IO = IGO)
io_top_names <- left_join(io_top, IO_names, by = "IO")

US <- state_mem %>% filter(state == 2) %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"),
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::scale_x_continuous(breaks = seq(min(year), max(year), by = 5))

China <- state_mem %>% filter(state == 710) %>% ggplot2::ggplot() +  #China
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::labs(y = "China Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

Taiwan <- state_mem %>% filter(state == 713) %>% ggplot2::ggplot() +  #Taiwan
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::labs(y = "Taiwan Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

Russia <- state_mem %>% filter(state == 365) %>% ggplot2::ggplot() +  #Russia
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::labs(y = "Russia Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))


South_Africa <- state_mem %>% filter(state == 560) %>% ggplot2::ggplot() +  #South Africa
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::labs(y = "South Africa Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

Vietnam <- state_mem %>% filter(state == 816) %>% ggplot2::ggplot() +  #Vietnam
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::labs(y = "Vietnam Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

Japan <- state_mem %>% filter(state == 740) %>% ggplot2::ggplot() +  #Japan
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) +
  ggplot2::labs(y = "Japan Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

Argentina <- state_mem %>% filter(state == 160) %>% ggplot2::ggplot() +  #Argentina
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::labs(y = "Argentina Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

pdf("state_mem.pdf", width = 13, height = 8)
gridExtra::grid.arrange(China, Taiwan, Russia, South_Africa, Japan, Argentina, ncol = 3)
dev.off()  # Close the PDF device

oecd <- io_mem %>% filter(IO == "OECD") %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "OECD Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

iaea <- io_mem %>% filter(IO == "IAEA") %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "IAEA Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

au <- io_mem %>% filter(IO == "AU") %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "AU Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

un <- io_mem %>% filter(IO == "UN") %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "UN Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

sco <- io_mem %>% filter(IO == "SCO") %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "SCO Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

wto <- io_mem %>% mutate(wto = ifelse(IO == "GATT" | IO == "WTO", 1, 0), overlap = ifelse(IO == "WTO" & year == 1994, 1, 0)) %>% 
  filter(wto == 1 & overlap == 0) %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "GATT/WTO Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))

pdf("io_mem.pdf", width = 13, height = 8)
gridExtra::grid.arrange(un, au, oecd, iaea, sco, wto, ncol = 3)
dev.off()  # Close the PDF device



oas <- io_mem %>% filter(IO == "OAS") %>% ggplot2::ggplot() + 
  ggplot2::geom_area(ggplot2::aes_string(y = "value", x = "year", fill="group"), 
                     stat="identity", position="stack") + 
  ggplot2::guides(fill=ggplot2::guide_legend(title="Group")) + 
  ggplot2::labs(y = "OAS Mixed-Membership") +
  ggplot2::scale_x_continuous(breaks = seq(min(state_mem$year), max(state_mem$year), by = 10))



