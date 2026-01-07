###################################
#### clean working environment ####
###################################

rm(list = ls())

###################
#### Load data ####
###################

load("~/data_wytham_tits.RData")

########################
#### Make Figure S1 ####
########################

ggplot(data = df_fecundity, aes(x=ClutchSize_observed - BroodSize_observed)) +
  geom_histogram()+
  xlab("Brood reduction (number of nestling)")+
  ylab("Number of broods")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))
