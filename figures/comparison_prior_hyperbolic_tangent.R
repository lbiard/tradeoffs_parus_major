###################################
#### clean working environment ####
###################################

rm(list = ls())

#######################
#### load packages ####
#######################

library(ggplot2)

###########################################################################
#### Compare prior distributions with hyperbolic tangent link function ####
###########################################################################

# "narrower" prior distributon, based on a Normal(0,0.5) 
narrow <- tanh(rnorm(100000, 0, 0.5))
# "wider" prior distributon, based on a Normal(0,1)
wide <- tanh(rnorm(100000, 0, 1))

data_plot <- as.data.frame(cbind(narrow, wide))

# Plot the density of both distribution on the "correlation" scale ([-1,1])

ggplot(data_plot, aes(x=narrow)) +
  geom_density(fill="red", alpha=0.5)+
  geom_density(aes(x=wide), fill="blue", alpha=0.5)+
  xlab("Value (i.e. correlation)")+
  ylab("Density")+
  coord_cartesian(clip = 'off')+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate("text", x = 0.4, y = 0.9, label = "N(0,0.5)", color="red")+
  annotate("text", x = -0.4, y = 0.9, label = "N(0,1)", color="blue")




