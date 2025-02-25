#######################
#### load packages ####
#######################

library(ggplot2)
library(tibble)
library(ggdist)
library(bayesplot)
library(patchwork)

########################################
#### Load dataset and  model output ####
########################################

dat_plot <- read.csv("~/My Drive/phd/chapter 3/wytham_wood/analysis/last_version/last_last_version/posterior_draws.txt", row.names=1, sep="")

df <- read.delim("~/My Drive/phd/chapter 3/wytham_wood/analysis/last_version/last_last_version/df.txt")
df_fecundity <- read.delim("~/My Drive/phd/chapter 3/wytham_wood/analysis/last_version/last_last_version/df_fecundity.txt")







####################################
####################################
#### Effects on trait variances ####
####################################
####################################

#####################################################
#### Among-individual variance in offspring mass ####
#####################################################

#### Density ####

x2.sim = seq(min(scale(df_fecundity$population_density)[,1]),
             max(scale(df_fecundity$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.1. + dat_plot$B_v.4.1. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.1.) + median(dat_plot$B_v.4.1.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_v1 <- function(x){(x*sd(df_fecundity$population_density))+mean(df_fecundity$population_density)}


plot_vof_density <- ggplot(plot.dat, aes(x = inv_fun_v1(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Population density")+
  ylab("Among mother variance in offspring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



#### temperature ####


x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.1. + dat_plot$B_v.2.1. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.1.) + median(dat_plot$B_v.2.1.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_v2 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature))+mean(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)}


plot_vof_temperature <- ggplot(plot.dat, aes(x = inv_fun_v2(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring temperature")+
  ylab("Among mother variance in offspring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))




#### precipitation ####


x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.1. + dat_plot$B_v.3.1. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.1.) + median(dat_plot$B_v.3.1.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_v3 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation))+mean(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)}


plot_vof_precipitation <- ggplot(plot.dat, aes(x = inv_fun_v3(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring precipitation")+
  ylab("Among mother variance in offspring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



#### beech mast index ####


x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.1. + dat_plot$B_v.5.1. * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.1.) + median(dat_plot$B_v.5.1.) * (x2.sim)) 
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)


plot_vof_beechmast <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Beech mast index")+
  ylab("Among mother variance in offspring mass")+
  scale_x_continuous(breaks=c(0,1,2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



combined_plot <- plot_vof_density+plot_vof_temperature+plot_vof_precipitation+plot_vof_beechmast & ylab(NULL) & ggtitle(NULL)

#ragg::agg_tiff("Figure S3.tiff", width = 8, height = 6, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
gridExtra::grid.arrange(gt, left = "Among mother variance in offspring mass", top = "Offpring mass variance")  
#dev.off()


##################################################
#### observation level variance in brood size ####
##################################################


#### Density ####

x2.sim = seq(min(scale(df_fecundity$population_density)[,1]),
             max(scale(df_fecundity$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.2. + dat_plot$B_v.4.2. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.2.) + median(dat_plot$B_v.4.2.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_vb1 <- function(x){(x*sd(df_fecundity$population_density))+mean(df_fecundity$population_density)}


plot_vb_density <- ggplot(plot.dat, aes(x = inv_fun_vb1(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Population density")+
  ylab("Observation level variance in brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



#### temperature ####


x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.2. + dat_plot$B_v.2.2. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.2.) + median(dat_plot$B_v.2.2.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_vb2 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature))+mean(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)}


plot_vb_temperature <- ggplot(plot.dat, aes(x = inv_fun_vb2(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring temperature")+
  ylab("Observation level variance in brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))




#### precipitation ####


x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.2. + dat_plot$B_v.3.2. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.2.) + median(dat_plot$B_v.3.2.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_vb3 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation))+mean(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)}


plot_vb_precipitation <- ggplot(plot.dat, aes(x = inv_fun_vb3(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring precipitation")+
  ylab("Observation level variance in brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



#### beech mast index ####


x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.2. + dat_plot$B_v.5.2. * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.2.) + median(dat_plot$B_v.5.2.) * (x2.sim)) 
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)


plot_vb_beechmast <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Beech mast index")+
  ylab("Observation level variance in brood size")+
  scale_x_continuous(breaks=c(0,1,2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



combined_plot <- plot_vb_density+plot_vb_temperature+plot_vb_precipitation+plot_vb_beechmast & ylab(NULL) & ggtitle(NULL)

#ragg::agg_tiff("Figure S4.tiff", width = 8, height = 6, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
gridExtra::grid.arrange(gt, left = "Observation level variance in brood size", top = "Brood size variance")  
#dev.off()




###################################################
#### observation level variance in recruitment ####
###################################################


#### Density ####

x2.sim = seq(min(scale(df_fecundity$population_density)[,1]),
             max(scale(df_fecundity$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.3. + dat_plot$B_v.4.3. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.3.) + median(dat_plot$B_v.4.3.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_vr1 <- function(x){(x*sd(df_fecundity$population_density))+mean(df_fecundity$population_density)}


plot_vr_density <- ggplot(plot.dat, aes(x = inv_fun_vr1(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Population density")+
  ylab("Observation level variance in recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



#### temperature ####


x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.3. + dat_plot$B_v.2.3. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.3.) + median(dat_plot$B_v.2.3.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_vr2 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature))+mean(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)}


plot_vr_temperature <- ggplot(plot.dat, aes(x = inv_fun_vr2(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring temperature")+
  ylab("Observation level variance in recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))




#### precipitation ####


x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.3. + dat_plot$B_v.3.3. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.3.) + median(dat_plot$B_v.3.3.) * (x2.sim))
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)
inv_fun_vr3 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation))+mean(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)}


plot_vr_precipitation <- ggplot(plot.dat, aes(x = inv_fun_vr3(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring precipitation")+
  ylab("Observation level variance in recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



#### beech mast index ####


x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_v.1.3. + dat_plot$B_v.5.3. * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_v.1.3.) + median(dat_plot$B_v.5.3.) * (x2.sim)) 
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)


plot_vr_beechmast <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Beech mast index")+
  ylab("Observation level variance in recruitment")+
  scale_x_continuous(breaks=c(0,1,2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



combined_plot <- plot_vr_density+plot_vr_temperature+plot_vr_precipitation+plot_vr_beechmast & ylab(NULL) & ggtitle(NULL)

#ragg::agg_tiff("Figure S5.tiff", width = 8, height = 6, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
gridExtra::grid.arrange(gt, left = "Observation level variance in recruitment", top = "Recruitment variance")  
#dev.off()



####################################
#### Figure posteriors variance ####
####################################


# change the "30000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 30 = 30000)

niter=60000 

df.posteriors <- data_frame(Submodel = c(rep("offspring mass variance", niter*4), rep("brood size variance", niter*4), rep("recruitment variance", niter*4))
                            , parameter = c(rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter)
                                            , rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter)
                                            , rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter))
                            , Posterior = c(dat_plot$B_v.4.1., dat_plot$B_v.2.1., dat_plot$B_v.3.1., dat_plot$B_v.5.1.
                                            , dat_plot$B_v.4.2., dat_plot$B_v.2.2., dat_plot$B_v.3.2., dat_plot$B_v.5.2.
                                            , dat_plot$B_v.4.3., dat_plot$B_v.2.3., dat_plot$B_v.3.3., dat_plot$B_v.5.3.))



df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("offspring mass variance", "brood size variance", "recruitment variance"))
df.posteriors$parameter <- factor(df.posteriors$parameter, levels = c("Beech mast index", "Spring precipitation", "Spring temperature", "Population density"))

submodels <- unique(df.posteriors$Submodel)
parameters <- unique(df.posteriors$parameter)

# get summaries of posteriors
df.summary <- expand.grid(Submodel = factor(x = submodels, levels = submodels, ordered = TRUE)
                          , parameter = factor(x = parameters, levels = parameters, ordered = TRUE)
                          , mean = as.numeric(NA)
                          , BCI50_upper = as.numeric(NA)
                          , BCI50_lower = as.numeric(NA)
                          , BCI89_upper = as.numeric(NA)
                          , BCI89_lower = as.numeric(NA)
                          , significant = as.logical(NA))


# get quantiles
for(i in 1:nrow(df.summary)){
  row <- which(df.posteriors$Submodel==as.character(df.summary$Submodel[i]) & df.posteriors$parameter==as.character(df.summary$parameter[i]))
  df.summary$mean[i] <- mean(df.posteriors$Posterior[row])
  df.summary$BCI50_upper[i] <- quantile(df.posteriors$Posterior[row], 0.75, na.rm = T)
  df.summary$BCI50_lower[i] <- quantile(df.posteriors$Posterior[row], 0.25, na.rm = T)
  df.summary$BCI89_upper[i] <- quantile(df.posteriors$Posterior[row], 0.945, na.rm = T)
  df.summary$BCI89_lower[i] <- quantile(df.posteriors$Posterior[row], 0.055, na.rm = T)
  df.summary$significant[i] <- ifelse(test = (df.summary$BCI89_lower[i]>0 & df.summary$BCI89_upper[i]>0) || (df.summary$BCI89_lower[i]<0 & df.summary$BCI89_upper[i]<0), yes = TRUE, no = FALSE)
}



# some plot settings, choose some pretty colors
dodge.width <- 0.7
colT1 <- "cornflowerblue"
colT2 <- "tomato3"
colT3 <- "olivedrab3"
colPlot <- "black"

# plot
p2 <-ggplot()+
  stat_halfeye(data = df.posteriors, aes(x = Posterior, y = parameter, fill = Submodel), color = NA, alpha = 0.2, position = position_dodge(width = dodge.width), normalize="xy", scale=0.55)+
  geom_point(data = df.summary, aes(x = mean, y = parameter, color = Submodel), size = 2, position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI89_lower, xmax = BCI89_upper, y = parameter, color = Submodel), size=0.6, linetype="solid", position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI50_lower, xmax = BCI50_upper, y = parameter, color = Submodel), size=1.5, linetype="solid", position = position_dodge(width = dodge.width))+
  scale_color_manual(values = c(colT1, colT2, colT3))+
  scale_fill_manual(values = c(colT1, colT2, colT3))+
  scale_alpha_manual(values = c(0.3,1))+
  geom_vline(xintercept = 0, linetype = "dashed", color = colPlot, size = 0.4) +
  scale_x_continuous(breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2))+
  coord_cartesian(xlim=c(-2.2,2.2), clip = "off")+
  ylab("")+
  theme_minimal()+
  theme(plot.margin = margin(20,5.5,5.5,5.5))+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )

#ragg::agg_tiff("Figure S6.tiff", width = 8, height = 7, units = "in", res = 300)
p2
#dev.off()




