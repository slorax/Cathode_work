#Data analysis for HEV program - Anh Vu

# Required packages to run the program ------------------------------------

library(tidyverse)
library(ggplot2)
library(data.table)
library(gridExtra)
library(cowplot)
library(ggpubr)

# Import maccor file and data clean up ------------------------------------

# Import excel file contains cell, laminate, and materials informations 

sample_list <- read.csv(file = "cell_info_03_06_2019.csv", stringsAsFactors = FALSE, header = TRUE)

# Function to import maccor files and add one column conatining full path to the file 

read_plus <- function(flnm) {
  fread(flnm, stringsAsFactors = FALSE)[ , 1:11] %>% 
    mutate(full_path = flnm)
}

# put all data in a dataframe and rename columns
data_join <- sample_list[ , "full_path"] %>% map_dfr(~ read_plus(.))
colnames(data_join) <- c("rec", 
                         "cycle", 
                         "Step", 
                         "test_time", 
                         "step_time", 
                         "capacity_Ah", 
                         "energy_Wh", 
                         "current_ampere", 
                         "potential", 
                         "state", 
                         "ES",
                         "full_path")
#adding electrode information to the echem results and calculate energy density, dqdv curve

HEV_data <- left_join(sample_list, data_join, by = "full_path") %>%
  mutate(capacity = capacity_Ah*1000*1000/active_mass_mg,
         energy = energy_Wh*1000000/active_mass_mg,
         current = current_ampere*1000,
         time = test_time/60,
         dqdV = c(0, diff(1000000*capacity_Ah/active_mass_mg)/diff(potential))) %>% 
  select(Cathode_ID, 
         precursor_notation, 
         Cell_ID, 
         Electrode_loading, 
         protocol, 
         Li_blending_ratio,
         Calcination_temp,
         Atmosphere,
         Maccor,
         cycle,
         Step,
         time,
         current,
         potential,
         capacity,
         energy,
         state,
         ES,
         dqdV)

#import data for cell tested using standard protocol
# calculate specific capacity of each cycle 

cell_data <- HEV_data %>%
  filter(state =="C"|state=="D") %>% 
  group_by(Cell_ID, precursor_notation, protocol, Li_blending_ratio, Calcination_temp, Atmosphere, cycle, state) %>% 
  slice(which.max(capacity)) %>% 
  select(Cell_ID, precursor_notation, protocol, Li_blending_ratio, Calcination_temp, Atmosphere, cycle, state, capacity, energy) %>% 
  mutate(avg_voltage = energy/capacity) %>% 
  arrange(Cell_ID, state)

#separate data for cells using different protocols
standard_protocols <- c("RNGC-v1", "RNGC-v2", "RNGC-v3", "1_RNGC_standard")
rate_protocols <- c("RNGC-RC-v1", "RNGC-RC-v2", "RNGC-RC-v3", "2_RNGC_rate")
camp_samples <- c("NMC 532", "NMC 622", "NMC 811", "NCA")

cell_standard <- cell_data %>% filter(protocol %in% standard_protocols)
cell_rate <- cell_data %>% filter(protocol %in% rate_protocols)

#calculate the hold discharge capacity for the first cycle

dchg_capacity_hold <-HEV_data %>% 
  filter((Step==6&ES==133)|(Step==7&ES==129), cycle==1) %>% 
  select(Cell_ID, Step, capacity) %>% 
  spread(key=Step, value = capacity) %>% 
  mutate(hold_dchg_capacity = .[[3]] - .[[2]]) %>% 
  select(Cell_ID, hold_dchg_capacity)

#calculate the hold chare capacity
chg_capacity_hold<-HEV_data %>% 
  filter((protocol %in% c("RNGC-v1", "RNGC-v2")&cycle==5&((Step==37&ES==133)|(Step==39&ES==129)))|
           ((protocol %in% c("RNGC-v3", "1_RNGC_standard")&cycle==7)&((Step==47&ES==133)|(Step==49&ES==129))), state=="C") %>% 
  select(Cell_ID, Step, capacity) %>%
  mutate( Step = factor(Step, levels = c(37, 39, 47, 49), labels = c("start hold", "end hold", "start hold", "end hold"))) %>% 
  spread(key = Step, value = capacity)%>% 
  mutate(hold_chg_capacity = .[[3]] - .[[2]]) %>% 
  select(Cell_ID, hold_chg_capacity)

#calculate polarization voltage for cell using old protocols
voltage_pol_old_protocol <- HEV_data %>% 
  filter((protocol %in% c("RNGC-v1", "RNGC-v2")&cycle==3), 
         (Step %in% c(15, 17, 19, 23, 25, 27)&ES==164)|
           (Step %in% c(16, 18, 20, 22, 24, 26, 28, 30)&ES==129)|
           (Step %in% c(21, 29)&ES==133)) %>%
  select(Cell_ID, Step, potential) %>% 
  spread(key = Step, value = potential) %>% 
  mutate(chg_pol_25 = (.[[2]] - .[[3]])*1000,
         chg_pol_50 = (.[[4]] - .[[5]])*1000,
         chg_pol_75 = (.[[6]] - .[[7]])*1000,
         chg_pol_100 = (.[[8]] - .[[9]])*1000,
         dchg_pol_25 = (.[[11]] - .[[10]])*1000,
         dchg_pol_50 = (.[[13]] - .[[12]])*1000,
         dchg_pol_75 = (.[[15]] - .[[14]])*1000,
         dchg_pol_100 = (.[[17]] - .[[16]])*1000) %>% 
  select(Cell_ID, chg_pol_25, chg_pol_50, chg_pol_75, chg_pol_100, dchg_pol_25, dchg_pol_50, dchg_pol_75, dchg_pol_100)

#calculate polarization voltage for cell using new protocols
voltage_pol_new_protocol <- HEV_data %>% 
  filter((protocol %in% c("RNGC-v3", "1_RNGC_standard")&cycle==4), 
         (Step %in% c(20, 22, 24, 28, 30, 32)&ES==164)|
           (Step %in% c(21, 23, 25, 27, 29, 31, 33)&ES==129)|(Step==26)&ES==133) %>%
  select(Cell_ID, Step, potential) %>% 
  spread(key = Step, value = potential) %>%
  mutate(chg_pol_25 = (.[[2]] - .[[3]])*1000,
         chg_pol_50 = (.[[4]] - .[[5]])*1000,
         chg_pol_75 = (.[[6]] - .[[7]])*1000,
         chg_pol_100 = (.[[8]] - .[[9]])*1000,
         dchg_pol_25 = (.[[11]] - .[[10]])*1000,
         dchg_pol_50 = (.[[13]] - .[[12]])*1000,
         dchg_pol_75 = (.[[15]] - .[[14]])*1000)%>% 
  select(Cell_ID, chg_pol_25, chg_pol_50, chg_pol_75, chg_pol_100, dchg_pol_25, dchg_pol_50, dchg_pol_75)

#joining two tables
voltage_pol <- full_join(voltage_pol_old_protocol, voltage_pol_new_protocol) %>% 
  arrange(desc(Cell_ID))

#final table contain hold chg, dchg capacity and polarization
hold_chg_dchg_pol <- sample_list %>% 
  filter(protocol %in% standard_protocols) %>% 
  full_join(dchg_capacity_hold) %>% 
  full_join(chg_capacity_hold) %>% 
  full_join(voltage_pol) %>% 
  select(-c(full_path, active_mass_mg, Active_loading, Electrode_loading, Maccor))

#transform the cell data frame to get capacity, energy in each column
#ungroup to remove "State" grouping variable (or change it to a fixed string, may not nescessary and depending on different R version)
cycle_standard_max <- max(cell_standard$cycle)

chg_cap_standard <- cell_standard %>% 
  filter(state=="C", cycle %in% c(1:cycle_standard_max)) %>%
  ungroup(state) %>% 
  mutate(cycle = factor(cycle, levels = c(1:cycle_standard_max), labels = paste("chg_cap_", 1:cycle_standard_max, sep = ""))) %>% 
  select(Cell_ID, cycle, capacity) %>%
  spread(key = cycle, value = capacity)


dchg_cap_standard <- cell_standard %>% 
  filter(state=="D", cycle %in% c(1:cycle_standard_max)) %>%
  ungroup(state) %>%
  mutate(cycle = factor(cycle, levels = c(1:cycle_standard_max), labels = paste("dchg_cap_", 1:cycle_standard_max, sep = ""))) %>% 
  select(Cell_ID, cycle, capacity) %>%
  spread(key = cycle, value = capacity) 
 
chg_ener_standard <- cell_standard %>% 
  filter(state=="C", cycle %in% c(1:cycle_standard_max)) %>%
  ungroup(state) %>%
  mutate(cycle = factor(cycle, levels = c(1:cycle_standard_max), labels = paste("chg_ener_", 1:cycle_standard_max, sep = ""))) %>% 
  select(Cell_ID, cycle, energy) %>%
  spread(key = cycle, value = energy) 

dchg_ener_standard <- cell_standard %>% 
  filter(state=="D", cycle %in% c(1:cycle_standard_max)) %>%
  ungroup(state) %>%
  mutate(cycle = factor(cycle, levels = c(1:cycle_standard_max), labels = paste("dchg_ener_", 1:cycle_standard_max, sep = ""))) %>% 
  select(Cell_ID, cycle, energy) %>%
  spread(key = cycle, value = energy) 

cell_standard_perfromance <- hold_chg_dchg_pol %>% 
  full_join(chg_cap_standard) %>% 
  full_join(dchg_cap_standard) %>% 
  full_join(chg_ener_standard) %>% 
  full_join(dchg_ener_standard) 
  

cell_standard_perfromance_avg <- cell_standard_perfromance %>% 
  filter(dchg_cap_1 > 50) %>%
  group_by(Cathode_ID, precursor_notation, protocol, Li_blending_ratio, Calcination_temp, Atmosphere) %>% 
  summarise_all(funs(mean, sd))

#remove temporary variables in the global environment
rm(chg_capacity_hold, 
   dchg_capacity_hold, 
   voltage_pol, 
   voltage_pol_new_protocol, 
   voltage_pol_old_protocol,
   chg_cap_standard,
   dchg_cap_standard,
   chg_ener_standard,
   dchg_ener_standard,
   hold_chg_dchg_pol,
   data_join)


# Data clean up for rate protocol 

chg_cap_rate <- cell_rate %>% 
  filter(state=="C", cycle %in% c(1:10)) %>%
  ungroup(state) %>% 
  mutate(cycle = factor(cycle, levels = c(1:10), labels = paste("chg_cap_", 1:10, sep = ""))) %>% 
  select(Cell_ID, cycle, capacity) %>%
  spread(key = cycle, value = capacity)


dchg_cap_rate <- cell_rate %>% 
  filter(state=="D", cycle %in% c(1:10)) %>%
  ungroup(state) %>%
  mutate(cycle = factor(cycle, levels = c(1:10), labels = paste("dchg_cap_", 1:10, sep = ""))) %>% 
  select(Cell_ID, cycle, capacity) %>%
  spread(key = cycle, value = capacity) 

chg_ener_rate <- cell_rate %>% 
  filter(state=="C", cycle %in% c(1:10)) %>%
  ungroup(state) %>%
  mutate(cycle = factor(cycle, levels = c(1:10), labels = paste("chg_ener_", 1:10, sep = ""))) %>% 
  select(Cell_ID, cycle, energy) %>%
  spread(key = cycle, value = energy) 

dchg_ener_rate <- cell_rate %>% 
  filter(state=="D", cycle %in% c(1:10)) %>%
  ungroup(state) %>%
  mutate(cycle = factor(cycle, levels = c(1:10), labels = paste("dchg_ener_", 1:10, sep = ""))) %>% 
  select(Cell_ID, cycle, energy) %>%
  spread(key = cycle, value = energy) 

cell_rate_perfromance <- sample_list %>% 
  filter(protocol %in% rate_protocols) %>% 
  full_join(chg_cap_rate) %>% 
  full_join(dchg_cap_rate) %>% 
  full_join(chg_ener_rate) %>% 
  full_join(dchg_ener_rate) %>% 
  select(-c("full_path"))


cell_rate_perfromance_avg <- cell_rate_perfromance %>% 
  filter(dchg_cap_1 > 50) %>%
  group_by(Cathode_ID, precursor_notation, protocol, Li_blending_ratio, Calcination_temp, Atmosphere) %>% 
  summarise_all(funs(mean, sd))

rm(chg_cap_rate,
   dchg_cap_rate,
   chg_ener_rate,
   dchg_ener_rate)

# calculate average for cells using rate protocol 


cell_rate_avg <- cell_rate %>% 
  ungroup (Cell_ID) %>% 
  select(- Cell_ID) %>%
  group_by(precursor_notation,
           protocol,
           Li_blending_ratio,
           Calcination_temp, 
           Atmosphere,
           cycle,
           state) %>% 
  summarise(cap_avg = mean(capacity),
            ener_avg = mean(energy),
            vol_avg = mean(avg_voltage),
            cap_sd = sd(capacity),
            ener_sd = sd(energy),
            vol_sd = sd(avg_voltage))

write.csv(cell_standard_perfromance, file = "cell_standard_performance.csv")
# Ploting section ---------------------------------------------------------

# standard protocol - effect of Li stoichimometry and atmosphere ----------------------
#change temp as needed
cell_standard_plotting <- cell_standard_perfromance %>% filter(dchg_cap_1 > 70, dchg_cap_6>70, Calcination_temp == 850)

chg_hold_plotting <- cell_standard_perfromance %>% filter((protocol %in% c("RNGC-v1", "RNGC-v2") & cycle == 5)|
                                                            (protocol %in% c("RNGC-v3", "1-RNGC-standard")&cycle == 7))
view(cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))

#capacity--------
g1<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_2, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 2.5V-4.2V', 
       x = 'Li/TM blending ratio', 
       y = 'Dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(65, 145))+
  theme_bw()

g2<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_6, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 2.5V-4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(75, 185), breaks = seq(75, 185, 10))+
  theme_bw()

g3<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_2, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 2.5V-4.2V ', 
       x = 'Li/TM blending ratio', 
       y = 'Dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(65, 145))+
  theme_bw()

g4<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_6, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 2.5V-4.5V ', 
       x = 'Li/TM blending ratio', 
       y = 'Dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(75, 185), breaks = seq(75, 185, 10))+
  theme_bw()

g5<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_2, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 2.5V-4.2V ', 
       x = 'Li/TM blending ratio', 
       y = 'Dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(65, 145))+
  theme_bw()

g6<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_6, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 2.5V-4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(75, 185), breaks = seq(75, 185, 10))+
  theme_bw()

ggarrange(g1, g3, g5, g2, g4, g6, 
                    ncol=3, 
                    nrow=2, 
                    common.legend = TRUE, 
                    legend = "bottom", 
                    labels = c("A", "B", "C", "D", "E", "F"))

rm(g1, g2, g3, g4, g5, g6)

#hold chg capacity and hold dchg capacity--------
g1<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_dchg_capacity, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(10, 45))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g2<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_chg_capacity, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 - 60h hold at 4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold chg capacity, mAh/g')+
  scale_y_continuous(limits = range(35, 90))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g3<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_dchg_capacity, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(10, 45))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g4<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_chg_capacity, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 - 60h hold at 4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold chg capacity, mAh/g')+
  scale_y_continuous(limits = range(35, 90))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g5<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_dchg_capacity, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(10, 45))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g6<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_chg_capacity, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 - 60h hold at 4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold chg capacity, mAh/g')+
  scale_y_continuous(limits = range(35, 90))+
  theme_bw()

figure <- ggarrange(g1, g3, g5, g2, g4, g6, 
                    ncol=3, 
                    nrow=2, 
                    common.legend = TRUE, 
                    legend = "bottom", 
                    labels = c("A", "B", "C", "D", "E", "F"))

annotate_figure(figure,
                top = text_grob("Effects of Li/TM ratios and atmosphere on hold chg/dchg capacity", color = "blue", face = "bold", size = 18))

rm(figure, g1, g2, g3, g4, g5, g6)

#kinetic dchg capacity--------
g1<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_dchg_capacity, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(10, 45))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g2<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_1, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold chg capacity, mAh/g')+
  scale_y_continuous(limits = range(90, 160))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g3<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_dchg_capacity, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(10, 45))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g4<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_1, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold chg capacity, mAh/g')+
  scale_y_continuous(limits = range(90, 160))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g5<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = hold_dchg_capacity, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold dchg capacity, mAh/g')+
  scale_y_continuous(limits = range(10, 45))+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g6<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_cap_1, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 - 20h hold at 2.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Hold chg capacity, mAh/g')+
  scale_y_continuous(limits = range(90, 160))+
  theme_bw()

figure <- ggarrange(g1, g3, g5, g2, g4, g6, 
                    ncol=3, 
                    nrow=2, 
                    common.legend = TRUE, 
                    legend = "bottom", 
                    labels = c("A", "B", "C", "D", "E", "F"))

annotate_figure(figure,
                top = text_grob("Effects of Li/TM ratios and atmosphere on kinetic dchg capaciity", color = "blue", face = "bold", size = 18))

rm(figure, g1, g2, g3, g4, g5, g6)


#polarization during chg and dchg-------------

g1<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = chg_pol_50, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 2.5V-4.2V', 
       x = 'Li/TM blending ratio', 
       y = 'Polarization, mV')+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g2<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_pol_50, fill=Atmosphere))+
  labs(title = 'Mn50Ni50 2.5V-4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Polarization, mV')+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g3<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = chg_pol_50, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 2.5V-4.2V', 
       x = 'Li/TM blending ratio', 
       y = 'Polarization, mV')+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g4<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_pol_50, fill=Atmosphere))+
  labs(title = 'Mn44Ni56 2.5V-4.5V', 
       x = 'Li/TM blending ratio', 
       y = 'Polarization, mV')+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g5<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = chg_pol_50, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 2.5V-4.2V', 
       x = 'Li/TM blending ratio', 
       y = 'Polarization, mV')+
  scale_fill_discrete(guide=FALSE)+
  theme_bw()

g6<-ggplot(data = cell_standard_plotting %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_boxplot(aes(x = factor(round(Li_blending_ratio, digits = 3)), y = dchg_pol_50, fill=Atmosphere))+
  labs(title = 'Mn41Ni54Co05 2.5V-4.2V', 
       x = 'Li/TM blending ratio', 
       y = 'Polarization, mV')+
  theme_bw()

figure <- ggarrange(g1, g3, g5, g2, g4, g6, 
          ncol=3, 
          nrow=2, 
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("A", "B", "C", "D", "E", "F"))

annotate_figure(figure,
                top = text_grob("Polarization of electrode at 50% SOD", color = "blue", face = "bold", size = 18))

rm(figure, g1, g2, g3, g4, g5, g6)

# Compare best samples from Mn50Ni50 and Mn41Ni54Co05 and NMC622 ----------

# select samples to be plotted --------------------------------------------
standard_samples_best <- cell_standard_perfromance %>% 
  group_by(precursor_notation, Atmosphere) %>% 
  slice(which.max(dchg_cap_6)) %>% 
  ungroup() %>% 
  select(Cell_ID, precursor_notation, Atmosphere, Calcination_temp, dchg_cap_6)
c("DO-HEV-139-08", "DO-HEV-148-08", "DO-HEV-142-15", "DO-HEV-146-14", "DO-HEV-622-01")
sample_standard_comparision <- c("DO-HEV-622-04", "AV-HEV-033-2", "AV-HEV-033-7")
#chg-dchg profiles 4.2V cut-off

ggplot() +
  geom_line(data = HEV_data %>% 
              filter(Cell_ID %in% sample_standard_comparision, cycle ==1, state=="C"), 
            aes(y = potential, 
                x = capacity, 
                color = precursor_notation,
                linetype = Atmosphere), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==1, state=="D"), 
            aes(y = potential, 
                x = capacity, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  labs(title = 'chg/dchg of Mn50Ni50 - Mn41Ni54Co05 - NMC 622 at 4.2V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'potential, V')+
  scale_y_continuous(limits = range(2.5, 4.25))+
  scale_x_continuous(limits = range(0, 200))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_colour_discrete(name="Composition")

ggplot() +
  geom_line(data = HEV_data %>% 
              filter(Cell_ID %in% sample_standard_comparision, cycle ==2, state=="C"), 
            aes(y = potential, 
                x = capacity, 
                color = precursor_notation,
                linetype = Atmosphere), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==2, state=="D"), 
            aes(y = potential, 
                x = capacity, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  labs(title = 'chg/dchg of Mn50Ni50 - Mn41Ni54Co05 - NMC 622 at 4.2V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'potential, V')+
  scale_y_continuous(limits = range(2.5, 4.25))+
  scale_x_continuous(limits = range(0, 200))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_colour_discrete(name="Composition")

#chg-dchg profiles 5th cycle
ggplot() +
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==6, state=="C"), 
            aes(y = potential, 
                x = capacity, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==6, state=="D"), 
            aes(y = potential, 
                x = capacity, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  labs(title = 'cchg/dchg of Mn50Ni50 - Mn41Ni54Co05 - NMC 622 at 4.5V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'potential, V')+
  scale_y_continuous(limits = range(2.5, 4.6))+
  scale_x_continuous(limits = range(0, 220))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_colour_discrete(name="Composition")


#dq/dV 1st cycle
ggplot() +
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==1, state=="C"), 
            aes(x = potential, 
                y = dqdV, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==1, state=="D"), 
            aes(x = potential, 
                y = dqdV, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  labs(title = 'dQ/dV of Mn50Ni50 - Mn41Ni54Co05 - NMC 622-4.2V cut-off', 
       x = 'potential, V', 
       y = 'dQ/dV, mAh/gV')+
  scale_y_continuous(limits = range(-700, 1800))+
  scale_x_continuous(limits = range(3.3, 4.25))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_colour_discrete(name="Composition")

ggplot() +
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==5, state=="C"), 
            aes(x = potential, 
                y = dqdV, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID %in% sample_standard_comparision, cycle ==5, state=="D"), 
            aes(x = potential, 
                y = dqdV, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  labs(title = 'dQ/dV of Mn50Ni50 - Mn41Ni54Co05 - NMC 622-4.5V cut-off', 
       x = 'potential, V', 
       y = 'dQ/dV, Ah/V')+
  scale_y_continuous(limits = range(-700, 1200))+
  scale_x_continuous(limits = range(3.3, 4.25))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_colour_discrete(name="Composition")

#dq/dV 6th cycle
ggplot() +
  geom_line(data = HEV_data %>% filter(Cell_ID %in% best_standard_samples, cycle ==5, state=="C", dqdV>0), 
            aes(x = potential, 
                y = dqdV, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID %in% best_standard_samples, cycle ==5, state=="D", dqdV<0), 
            aes(x = potential, 
                y = dqdV, 
                color = precursor_notation,
                linetype=Atmosphere), lwd = 1)+
  labs(title = 'dQ/dV of Mn50Ni50 - Mn41Ni54Co05 - NMC 622-4.5V cut-off', 
       x = 'potential, V', 
       y = 'dQ/dV, mAh/gV')+
  scale_y_continuous(limits = range(-600, 800))+
  scale_x_continuous(limits = range(3.3, 4.5))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_colour_discrete(name="Composition")
#compare ANL samples to NMC 622 from camp----------
cathode_composition <-c("Mn50Ni50", "Mn44Ni56", "Mn41Ni54Co05", "NMC 622")

# 4.2V cut-off -------------
#1st cycle with 20h hold at 2.5V
ggplot(data = cell_standard_perfromance_avg %>% filter(precursor_notation %in% cathode_composition)) +
  geom_point(aes(x = hold_dchg_capacity_mean, y = dchg_cap_1_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  labs(title = '20h discharge-hold at 2.5V', 
       x = 'Discharge-hold capacity, mAh/g', 
       y = 'Discharge capacity, mAh/g')+
  scale_y_continuous(limits = range(75, 200))+
  scale_x_continuous(limits = range(10, 45))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")

#2nd cycle without hold
ggplot(data = cell_standard_perfromance_avg %>% filter(precursor_notation %in% cathode_composition)) +
  geom_jitter(aes(y = dchg_cap_2_mean, x = dchg_cap_2_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  labs(title = '2.5V-4.2V - no holding step', 
       x = 'Charge capacity, mAh/g', 
       y = 'Discharge capacity, mAh/g')+
  scale_y_continuous(limits = range(50, 180))+
  scale_x_continuous(limits = range(50, 180))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")

#4.5V cut-off ---------

#60h hold at 4.5V
ggplot() +
  geom_point(data = cell_standard_perfromance_avg %>% 
               filter(precursor_notation %in% cathode_composition, protocol %in% c("RNGC-v1", "RNGC-v2")), 
             aes(x = hold_chg_capacity_mean, y = chg_cap_5_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  geom_point(data = cell_standard_perfromance_avg %>% 
               filter(precursor_notation %in% cathode_composition, protocol %in% c("RNGC-v3", "1_RNGC_standard")), 
             aes(x = hold_chg_capacity_mean, y = chg_cap_7_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  labs(title = '60h hold at 4.5V', 
       x = 'charge-hold capacity, mAh/g', 
       y = 'Charge capacity, mAh/g')+
  scale_y_continuous(limits = range(125, 225))+
  scale_x_continuous(limits = range(0, 80))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")

#2.5-4.5V cut-off 6th cycle
ggplot() +
  geom_point(data = cell_standard_perfromance_avg %>% 
               filter(precursor_notation %in% cathode_composition, protocol %in% c("RNGC-v1", "RNGC-v2")), 
             aes(x = dchg_cap_4_mean, y = chg_cap_4_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  geom_point(data = cell_standard_perfromance_avg %>% 
               filter(precursor_notation %in% cathode_composition, protocol %in% c("RNGC-v3", "1_RNGC_standard")), 
             aes(x = dchg_cap_6_mean, y = chg_cap_6_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  labs(title = '4.5V-cut off without 60h hold at 4.5V', 
       x = 'discharge capacity, mAh/g', 
       y = 'charge capacity, mAh/g')+
  scale_y_continuous(limits = range(125, 225))+
  scale_x_continuous(limits = range(75, 225))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")

#2.5-4.5V cut-off 5th cycle
ggplot() +
  geom_point(data = cell_standard_perfromance_avg %>% 
               filter(precursor_notation %in% cathode_composition, protocol %in% c("RNGC-v1", "RNGC-v2")), 
             aes(x = dchg_cap_4_mean, y = chg_cap_4_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  geom_point(data = cell_standard_perfromance_avg %>% 
               filter(precursor_notation %in% cathode_composition, protocol %in% c("RNGC-v3", "1_RNGC_standard")), 
             aes(x = dchg_cap_5_mean, y = chg_cap_5_mean, color=precursor_notation, shape = Atmosphere), size=3)+
  labs(title = '4.5V - cutoff', 
       x = 'Discharge capacity, mAh/g', 
       y = 'Charge capacity, mAh/g')+
  scale_y_continuous(limits = range(75, 225))+
  scale_x_continuous(limits = range(75, 225))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")

# plots for samples tested with rate protocol -----------------------------
sample_rate_air <- cell_rate_avg %>% filter(Calcination_temp==850, Atmosphere =="Air", Li_blending_ratio > 0.955)
sample_rate_oxygen <- cell_rate_avg %>% filter(Calcination_temp==850, Atmosphere =="Oxygen", Li_blending_ratio > 0.955)

# #Rate - effects of composition ------------------------------------------

g1<-ggplot(data = sample_rate_air %>% filter(round(Li_blending_ratio, digits = 3) ==1.015))+
  geom_point(aes(x = cycle, y = cap_avg, color=precursor_notation))+
  labs(title = "Rate performance of samples calcined under air",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_discrete(name="Composition")+
  scale_color_brewer(palette="Dark2")


g2<-ggplot(data = sample_rate_oxygen %>% filter(round(Li_blending_ratio, digits = 3) ==1.015))+
  geom_point(aes(x = cycle, y = cap_avg, color=precursor_notation))+
  labs(title = "Rate performance of sample calcined under oxygen",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Composition")

grid.arrange(g1, g2, ncol=1, nrow=2)

rm(g1, g2)

#Rate - effect of Li stoichiometry (change temperature, atmosphe --------

g1<-ggplot(data = sample_rate_air %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_point(aes(x = cycle, y = cap_avg, color=factor(round(Li_blending_ratio, digits = 3))))+
  labs(title = "Rate performance of Mn50Ni50 calcined under air",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Li/TM ratio")

g2<-ggplot(data = sample_rate_air %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_point(aes(x = cycle, y = cap_avg, color=factor(round(Li_blending_ratio, digits = 3))))+
  labs(title = "Rate performance of Mn44Ni56 calcined under air",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Li/TM ratio")

g3<-ggplot(data = sample_rate_air %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_point(aes(x = cycle, y = cap_avg, color=factor(round(Li_blending_ratio, digits = 3))))+
  labs(title = "Rate performance of Mn41Ni54Co05 calcined under air",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Li/TM ratio")

g4<-ggplot(data = sample_rate_oxygen %>% filter(precursor_notation=="Mn50Ni50"))+
  geom_point(aes(x = cycle, y = cap_avg, color=factor(round(Li_blending_ratio, digits = 3))))+
  labs(title = "Rate performance of Mn50Ni50 calcined under Oxygen",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Li/TM ratio")

g5<-ggplot(data = sample_rate_oxygen %>% filter(precursor_notation=="Mn44Ni56"))+
  geom_point(aes(x = cycle, y = cap_avg, color=factor(round(Li_blending_ratio, digits = 3))))+
  labs(title = "Rate performance of Mn44Ni56 calcined under Oxygen",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Li/TM ratio")

g6<-ggplot(data = sample_rate_oxygen %>% filter(precursor_notation=="Mn41Ni54Co05"))+
  geom_point(aes(x = cycle, y = cap_avg, color=factor(round(Li_blending_ratio, digits = 3))))+
  labs(title = "Rate performance of Mn41Ni54Co05 calcined under Oxygen",
       x = " Cycle",
       y = "Capacity, mAh/g")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_brewer(palette="Dark2")+
  scale_color_discrete(name="Li/TM ratio")

grid.arrange(g1, g3, g5, g2, g4, g6, ncol=3, nrow=2)

rm(g1, g2, g3, g4, g5, g6)

calcined_list <- sample_list %>% 
  filter(!precursor_notation %in% camp_samples) %>% 
  group_by(precursor_notation,
           Li_blending_ratio,
           Atmosphere,
           Calcination_temp) %>% 
  summarise() %>% 
  arrange(desc(precursor_notation)) %>% write.csv(file = "Calcination list.csv")





#low cut-off voltage 
ggplot()+
  geom_jitter(data = cell_standard_avg %>% filter(cycle==2, !precursor_notation %in% camp_samples, Li_blending_ratio > 0.9), 
             aes(x=cap_avg, y = energy_avg, shape= Atmosphere, color = factor(round(Li_blending_ratio, digits = 3))), size=3)+
  labs(title = 'Capacity vs. energy at 4.2V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'Energy, Wh/g')+
  scale_y_continuous(limits = range(200, 700))+
  scale_x_continuous(limits = range(60, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  scale_colour_discrete(name="Li/TM blending ratio",
                        labels = c("0.985",
                                   "1.015",
                                   "1.025",
                                   "1.07",
                                   "1.12"))+
  scale_color_brewer(palette="Dark2")

ggplot()+
  geom_jitter(data = cell_standard_avg %>% filter(cycle==2, precursor_notation =="Mn41Ni54Co05", Atmosphere=="Air"), 
              aes(x=cap_avg, y = energy_avg, shape= Atmosphere, color = factor(round(Li_blending_ratio, digits = 3))), size=3)+
  labs(title = 'Capacity vs. energy at 4.2V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'Energy, Wh/g')+
  scale_y_continuous(limits = range(200, 700))+
  scale_x_continuous(limits = range(60, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  scale_colour_discrete(name="Li/TM blending ratio",
                        labels = c("0.985",
                                   "1.015",
                                   "1.025",
                                   "1.07",
                                   "1.12"))

ggplot()+
  geom_jitter(data = cell_standard_avg %>% filter(cycle==2), 
             aes(x=cap_avg, y = energy_avg, shape= Atmosphere, color = precursor_notation), size=3)+
  labs(title = 'Capacity vs. energy at 4.2V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'Energy, Wh/g')+
  scale_y_continuous(limits = range(200, 700))+
  scale_x_continuous(limits = range(60, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  scale_colour_discrete(name="Composition")
                        

#high cut-off voltage
ggplot()+
  geom_point(data = cell_standard_avg %>% filter(cycle==5, !precursor_notation %in% camp_samples), 
             aes(x=cap_avg, y = energy_avg, shape= Atmosphere, color = factor(round(Li_blending_ratio, digits = 3))), size=3)+
  labs(title = 'Capacity vs. energy of at 4.5V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'Energy, Wh/g')+
  scale_y_continuous(limits = range(450, 900))+
  scale_x_continuous(limits = range(120, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  scale_colour_discrete(name="Li/TM blending ratio",
                        labels = c("0.985",
                                   "1.015",
                                   "1.025",
                                   "1.07",
                                   "1.12"))

ggplot()+
  geom_point(data = cell_standard_avg %>% filter(cycle==5), 
             aes(x=cap_avg, y = energy_avg, shape= Atmosphere, color = precursor_notation), size=3)+
  labs(title = 'Capacity vs. energy at 4.5V cut-off', 
       x = 'Capacity, mAh/g', 
       y = 'Energy, Wh/g')+
  scale_y_continuous(limits = range(450, 900))+
  scale_x_continuous(limits = range(110, 225))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  scale_colour_discrete(name="Composition")



high_vol_cap <-  cell_summary %>% 
  filter((protocol=="RNGC-v1"&cycle==5)|(protocol=="RNGC-v2"&cycle==6)|(protocol=="RNGC-v3"&cycle==6), state=="D") %>% 
  arrange(desc(capacity))

#remove bad cells and take average for cells of the same samples tested with rate protocol

cell_rate_ids <- cell_rate %>% filter(cycle==2, capacity > 70, state=="D") 
cell_rate_good <- cell_rate %>% filter(Cell_ID %in% cell_rate_ids$Cell_ID)
cell_rate_avg <- cell_rate_good %>% 
  group_by(precursor_notation, protocol, Li_blending_ratio, Calcination_temp, Atmosphere, cycle, state) %>% 
  summarise(cap_avg = mean(capacity),
            cap_sd = sd(capacity),
            energy_avg = mean(energy),
            energy_sd = sd(energy),
            voltage_avg = mean(avg_voltage),
            voltage_sd = sd(avg_voltage))


ggplot(data = cell_standard_avg %>% filter(cap_sd < 10))+
  geom_bar(aes(x = cap_sd), na.rm = TRUE, binwidth = 0.5)


# chg-dchg profiles -------------------------------------------------------

#single sample multiple cycles

ggplot() +
  geom_line(data = HEV_data %>% filter(Cell_ID =="Mn41Ni54Co05", cycle %in% cycle_plotting, state=="C") %>% group_by(cycle), 
            aes(y = potential, x = capacity, color=factor(cycle)), lwd = 1)+
  geom_line(data = HEV_data %>% filter(Cell_ID =="Mn41Ni54Co05", cycle %in% cycle_plotting, state=="D") %>% group_by(cycle), 
            aes(y = potential, x = capacity, color=factor(cycle)), lwd = 1)+
  labs(title = 'chg/dchg of 1st cycle', 
       x = 'Capacity, mAh/g', 
       y = 'potential, V')+
  scale_y_continuous(limits = range(3, 4.5))+
  scale_x_continuous(limits = range(0, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()

#single cycle multiple samples

ggplot() +
  geom_line(data = data_plotting %>% filter(cycle ==3, state=="C"), 
            aes(y = potential, x = capacity, color = Cell_ID), lwd = 1)+
  geom_line(data = data_plotting %>% filter(cycle ==3
                                            , state=="D"), 
            aes(y = potential, x = capacity, color = Cell_ID), lwd = 1)+
  labs(title = 'chg/dchg of 1st cycle', 
       x = 'Capacity, mAh/g', 
       y = 'potential, V')+
  scale_y_continuous(limits = range(3, 4.5))+
  scale_x_continuous(limits = range(0, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()

#multiple cycles of different samples put into grid

ggplot() +
  geom_line(data = data_plotting %>% filter(state=="C"), 
            aes(y = potential, x = capacity, color = Cell_ID), lwd = 1)+
  geom_line(data = data_plotting %>% filter(state=="D"), 
            aes(y = potential, x = capacity, color = Cell_ID), lwd = 1)+
  labs(title = 'chg/dchg of 1st cycle', 
       x = 'Capacity, mAh/g', 
       y = 'potential, V')+
  scale_y_continuous(limits = range(2.5, 4.5))+
  scale_x_continuous(limits = range(0, 200))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  facet_grid(rows = vars(cycle))


# dqdv plots --------------------------------------------------------------

#single cycle multiple samples

ggplot() +
  geom_line(data = data_plotting %>% filter(cycle ==1, state=="C"), 
            aes(x = potential, y = dqdV, color = Cell_ID), lwd = 1)+
  geom_line(data = data_plotting %>% filter(cycle ==1, state=="D"), 
            aes(x = potential, y = dqdV, color = Cell_ID), lwd = 1)+
  labs(title = 'dq/dV of 1st cycle', 
       x = 'potential, V', 
       y = 'dq/dV')+
  scale_y_continuous(limits = range(-0.01, 0.025))+
  scale_x_continuous(limits = range(3.25, 4.25))+
  theme(text = element_text(size = 16))+
  theme_bw()

#multiple cycles of different samples put into grid
ggplot() +
  geom_line(data = data_plotting %>% filter(state=="C"), 
            aes(x = potential, y = dqdV, color = Cell_ID), lwd = 1)+
  geom_line(data = data_plotting %>% filter(state=="D"), 
            aes(x = potential, y = dqdV, color = Cell_ID), lwd = 1)+
  labs(title = 'dq/dV of multiple cycles', 
       x = 'potential, V', 
       y = 'dq/dV')+
  scale_y_continuous(limits = range(-0.01, 0.025))+
  scale_x_continuous(limits = range(3.25, 4.25))+
  theme(text = element_text(size = 16))+
  theme_bw()+
  facet_grid(cols = vars(cycle))

ggplot(data=HEV_data %>% filter(Cell_ID=="AV-HEV-023-1"))+
  geom_line(aes(x = time, y = potential))
cycle4=HEV_data %>% filter(Cell_ID=="DO-HEV-622-01", cycle==4) %>% 
  select(state, ES, potential, cycle, Step)

short_table<-cell_standard_perfromance %>% 
  filter(!precursor_notation %in% camp_samples) %>%  
  select(Cell_ID,
         precursor_notation, 
         Calcination_temp, 
         Atmosphere, 
         Li_blending_ratio, 
         dchg_cap_1,
         dchg_cap_2,
         dchg_cap_5, 
         dchg_cap_6, 
         dchg_cap_7,
         hold_dchg_capacity,
         hold_chg_capacity) %>% 
  arrange(desc(dchg_cap_6))

sample_list %>% filter(protocol %in% standard_protocols) %>% 
  select(precursor_notation, Calcination_temp, Atmosphere, Li_blending_ratio) %>% write.csv(file = "cell_list_1.csv") 

