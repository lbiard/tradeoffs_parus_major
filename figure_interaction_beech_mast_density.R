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

dat_plot <- read.csv("~/posterior_draws_interaction.txt", row.names=1, sep="")

df <- read.delim("~/df.txt")
df_fecundity <- read.delim("~/df_fecundity.txt")

###########################
###########################
#### Figure posteriors ####
###########################
###########################

# change the "60000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 60 = 60000)

niter=60000 

df.posteriors <- data_frame(Submodel = c(rep("offspring mass - brood size", niter*6), rep("offpsring mass - recruitment", niter*6), rep("brood size - recruitment", niter*6))
                            , parameter = c(rep("Intercept", niter), rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Density * beech mast", niter)
                                            , rep("Intercept", niter), rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Density * beech mast", niter)
                                            , rep("Intercept", niter), rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Density * beech mast", niter))
                            , Posterior = c(dat_plot$B_cpc.1.1., dat_plot$B_cpc.4.1., dat_plot$B_cpc.2.1., dat_plot$B_cpc.3.1., dat_plot$B_cpc.5.1., dat_plot$B_cpc.6.1.
                                            , dat_plot$B_cpc.1.2., dat_plot$B_cpc.4.2., dat_plot$B_cpc.2.2., dat_plot$B_cpc.3.2., dat_plot$B_cpc.5.2., dat_plot$B_cpc.6.2.
                                            , dat_plot$B_cpc.1.3., dat_plot$B_cpc.4.3., dat_plot$B_cpc.2.3., dat_plot$B_cpc.3.3., dat_plot$B_cpc.5.3., dat_plot$B_cpc.6.3.))


df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("offspring mass - brood size", "offpsring mass - recruitment", "brood size - recruitment"))
df.posteriors$parameter <- factor(df.posteriors$parameter, levels = c("Density * beech mast", "Beech mast index", "Spring precipitation", "Spring temperature", "Population density", "Intercept"))

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
colCov <- "olivedrab3"
colPlot <- "black"

# plot
forest_plot <-ggplot()+
  stat_halfeye(data = df.posteriors, aes(x = Posterior, y = parameter, fill = Submodel), color = NA, alpha = 0.2, position = position_dodge(width = dodge.width), normalize="xy", scale=0.55)+
  geom_point(data = df.summary, aes(x = mean, y = parameter, color = Submodel), size = 2, position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI89_lower, xmax = BCI89_upper, y = parameter, color = Submodel), size=0.6, linetype="solid", position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI50_lower, xmax = BCI50_upper, y = parameter, color = Submodel), size=1.5, linetype="solid", position = position_dodge(width = dodge.width))+
  scale_color_manual(values = c(colT1, colT2, colCov))+
  scale_fill_manual(values = c(colT1, colT2, colCov))+
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

#ragg::agg_tiff("Figure S8.tiff", width = 8, height = 7, units = "in", res = 300)
forest_plot
#dev.off()



#######################################
#######################################
#### Effect on traits correlations ####
#######################################
#######################################

#####################################
#### offspring mass - brood size ####
#####################################

#### density ####

# sequence from scaled minimum to maximum value in x 
x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
             max(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
             by =  0.2) 

# low beech mast (0)
int.sim.low <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.low[, i] <- tanh(dat_plot$B_cpc.1.1. + dat_plot$B_cpc.4.1. * (x2.sim[i]) + dat_plot$B_cpc.6.1. * (x2.sim[i]) * 0) 
}

# high beech mast (2)
int.sim.high <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.high[, i] <- tanh(dat_plot$B_cpc.1.1. + dat_plot$B_cpc.4.1. * (x2.sim[i]) + dat_plot$B_cpc.6.1. * (x2.sim[i]) * 2) 
}


# calculate quantiles of predictions (low beech mast)
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.1.) + median(dat_plot$B_cpc.4.1.) * (x2.sim) + median(dat_plot$B_cpc.6.1.) * (x2.sim) * 0) 
bayes.c.eff.lower <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.low <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# calculate quantiles of predictions (high beech mast)
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.1.) + median(dat_plot$B_cpc.4.1.) * (x2.sim) + median(dat_plot$B_cpc.6.1.) * (x2.sim) * 2) 
bayes.c.eff.lower <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.high <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

inv_fun_p <- function(x){(x*sd(unique(df[,c("BreedingSeason","population_density")])$population_density))+mean(unique(df[,c("BreedingSeason","population_density")])$population_density)}


p <- ggplot(plot.dat.low, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#009E73", alpha = 0.1)+
  geom_line(color = "#009E73", size = 1.8, alpha=0.6)+
  
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#E69F00", alpha = 0.1)+
  geom_line(data=plot.dat.high, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean), color = "#E69F00", size = 1.8, alpha=0.6)+
  
  ylim(-1,1)+
  ggtitle("Offpsring mass - Brood size")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

p

# annotation
mesure_pol <- data.frame(x1 = c(-2.7,0.3),
                         x2 = c(-2.3,0.7),
                         cat = c(1:2),
                         catnames = c("Low beech mast","High beech mast"))

ann <- ggplot(data = mesure_pol) +
  geom_rect(aes(xmin = x1,
                xmax = x2,
                ymin = 0.095,
                ymax = 0.105,
                fill = as.factor(cat)),
            fill = c("#009E73", "#E69F00"),
            color = "black",
            size = 0.3) +
  geom_text(aes(x = x2+1.15, y = 0.1, label = catnames),
            #vjust = .8, 
            fontface = "bold", color = "black") +
  coord_cartesian(xlim = c(-2.7085, 2.861)
                  , ylim = c(0.065, 0.14)) +
  theme_void()
ann


#### beech mast ####

# sequence from minimum to maximum value in x 
x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

# low density (-1SD)
int.sim.low <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.low[, i] <- tanh(dat_plot$B_cpc.1.1. + dat_plot$B_cpc.5.1. * (x2.sim[i]) + dat_plot$B_cpc.6.1. * (x2.sim[i]) * -1) 
}

# high density (+1SD)
int.sim.high <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.high[, i] <- tanh(dat_plot$B_cpc.1.1. + dat_plot$B_cpc.5.1. * (x2.sim[i]) + dat_plot$B_cpc.6.1. * (x2.sim[i]) * 1) 
}


# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.1.) + median(dat_plot$B_cpc.5.1.) * (x2.sim) + median(dat_plot$B_cpc.6.1.) * (x2.sim) * -1)
bayes.c.eff.lower <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.low <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.1.) + median(dat_plot$B_cpc.5.1.) * (x2.sim) + median(dat_plot$B_cpc.6.1.) * (x2.sim) * 1)
bayes.c.eff.lower <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.high <- data.frame(x2.sim, bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

q <- ggplot(plot.dat.low, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#009E73", alpha = 0.1)+
  geom_line(color = "#009E73", size = 1.8, alpha=0.6)+
  
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#E69F00", alpha = 0.1)+
  geom_line(data=plot.dat.high, aes(x = x2.sim, y = bayes.c.eff.mean), color = "#E69F00", size = 1.8, alpha=0.6)+
  
  ylim(-1,1)+
  ggtitle("Offpsring mass - Brood size")+
  xlab("Beech mast index")+
  ylab("Observation-level correlation")+
  scale_x_continuous(breaks=c(0,1,2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
q


# annotation
mesure_pol_2 <- data.frame(x1 = c(-2.7,0.3),
                         x2 = c(-2.3,0.7),
                         cat = c(1:2),
                         catnames = c("Low density","High density"))

ann_2 <- ggplot(data = mesure_pol_2) +
  geom_rect(aes(xmin = x1,
                xmax = x2,
                ymin = 0.095,
                ymax = 0.105,
                fill = as.factor(cat)),
            fill = c("#009E73", "#E69F00"),
            color = "black",
            size = 0.3) +
  geom_text(aes(x = x2+1.15, y = 0.1, label = catnames),
            #vjust = .8, 
            fontface = "bold", color = "black") +
  coord_cartesian(xlim = c(-2.7085, 2.861)
                  , ylim = c(0.065, 0.14)) +
  theme_void()
ann_2


#####################################
#### ofspring mass - recruitment ####
#####################################

#### density ####

# sequence from scaled minimum to maximum value in x 
x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
             max(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
             by =  0.2) 

# low beech mast (0)
int.sim.low <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.low[, i] <- tanh(dat_plot$B_cpc.1.2. + dat_plot$B_cpc.4.2. * (x2.sim[i]) + dat_plot$B_cpc.6.2. * (x2.sim[i]) * 0) 
}

# high beech mast (2)
int.sim.high <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.high[, i] <- tanh(dat_plot$B_cpc.1.2. + dat_plot$B_cpc.4.2. * (x2.sim[i]) + dat_plot$B_cpc.6.2. * (x2.sim[i]) * 2) 
}


# calculate quantiles of predictions (low beech mast)
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.2.) + median(dat_plot$B_cpc.4.2.) * (x2.sim) + median(dat_plot$B_cpc.6.2.) * (x2.sim) * 0) 
bayes.c.eff.lower <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.low <- data.frame(x2.sim, bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# calculate quantiles of predictions (high beech mast)
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.2.) + median(dat_plot$B_cpc.4.2.) * (x2.sim) + median(dat_plot$B_cpc.6.2.) * (x2.sim) * 2) 
bayes.c.eff.lower <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.high <- data.frame(x2.sim, bayes.c.eff.mean,
                            bayes.c.eff.lower, bayes.c.eff.upper,
                            bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                            bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                            bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                            bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

inv_fun_p <- function(x){(x*sd(unique(df[,c("BreedingSeason","population_density")])$population_density))+mean(unique(df[,c("BreedingSeason","population_density")])$population_density)}


p2 <- ggplot(plot.dat.low, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#009E73", alpha = 0.1)+
  geom_line(color = "#009E73", size = 1.8, alpha=0.6)+
  
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#E69F00", alpha = 0.1)+
  geom_line(data=plot.dat.high, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean), color = "#E69F00", size = 1.8, alpha=0.6)+
  
  ylim(-1,1)+
  ggtitle("Offpsring mass - Recruitment")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

p2


#### beech mast ####

# sequence from minimum to maximum value in x 
x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

# low density (-1SD)
int.sim.low <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.low[, i] <- tanh(dat_plot$B_cpc.1.2. + dat_plot$B_cpc.5.2. * (x2.sim[i]) + dat_plot$B_cpc.6.2. * (x2.sim[i]) * -1) 
}

# high density (+1SD)
int.sim.high <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.high[, i] <- tanh(dat_plot$B_cpc.1.2. + dat_plot$B_cpc.5.2. * (x2.sim[i]) + dat_plot$B_cpc.6.2. * (x2.sim[i]) * 1) 
}


# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.2.) + median(dat_plot$B_cpc.5.2.) * (x2.sim) + median(dat_plot$B_cpc.6.2.) * (x2.sim) * -1)
bayes.c.eff.lower <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.low <- data.frame(x2.sim, bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.2.) + median(dat_plot$B_cpc.5.2.) * (x2.sim) + median(dat_plot$B_cpc.6.2.) * (x2.sim) * 1)
bayes.c.eff.lower <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.high <- data.frame(x2.sim, bayes.c.eff.mean,
                            bayes.c.eff.lower, bayes.c.eff.upper,
                            bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                            bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                            bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                            bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

q2 <- ggplot(plot.dat.low, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#009E73", alpha = 0.1)+
  geom_line(color = "#009E73", size = 1.8, alpha=0.6)+
  
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#E69F00", alpha = 0.1)+
  geom_line(data=plot.dat.high, aes(x = x2.sim, y = bayes.c.eff.mean), color = "#E69F00", size = 1.8, alpha=0.6)+
  
  ylim(-1,1)+
  ggtitle("Offpsring mass - Recruitment")+
  xlab("Beech mast index")+
  ylab("Observation-level correlation")+
  scale_x_continuous(breaks=c(0,1,2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
q2


##################################
#### brood size - recruitment ####
##################################

#### density ####

# sequence from scaled minimum to maximum value in x 
x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
             max(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
             by =  0.2) 

# low beech mast (0)
int.sim.low <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.low[, i] <- tanh(dat_plot$B_cpc.1.3. + dat_plot$B_cpc.4.3. * (x2.sim[i]) + dat_plot$B_cpc.6.3. * (x2.sim[i]) * 0) 
}

# high beech mast (2)
int.sim.high <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.high[, i] <- tanh(dat_plot$B_cpc.1.3. + dat_plot$B_cpc.4.3. * (x2.sim[i]) + dat_plot$B_cpc.6.3. * (x2.sim[i]) * 2) 
}


# calculate quantiles of predictions (low beech mast)
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.3.) + median(dat_plot$B_cpc.4.3.) * (x2.sim) + median(dat_plot$B_cpc.6.3.) * (x2.sim) * 0) 
bayes.c.eff.lower <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.low <- data.frame(x2.sim, bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# calculate quantiles of predictions (high beech mast)
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.3.) + median(dat_plot$B_cpc.4.3.) * (x2.sim) + median(dat_plot$B_cpc.6.3.) * (x2.sim) * 2) 
bayes.c.eff.lower <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.high <- data.frame(x2.sim, bayes.c.eff.mean,
                            bayes.c.eff.lower, bayes.c.eff.upper,
                            bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                            bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                            bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                            bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

inv_fun_p <- function(x){(x*sd(unique(df[,c("BreedingSeason","population_density")])$population_density))+mean(unique(df[,c("BreedingSeason","population_density")])$population_density)}


p3 <- ggplot(plot.dat.low, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#009E73", alpha = 0.1)+
  geom_line(color = "#009E73", size = 1.8, alpha=0.6)+
  
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#E69F00", alpha = 0.1)+
  geom_line(data=plot.dat.high, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean), color = "#E69F00", size = 1.8, alpha=0.6)+
  
  ylim(-1,1)+
  ggtitle("Brood size - Recruitment")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

p3


#### beech mast ####

# sequence from minimum to maximum value in x 
x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

# low density (-1SD)
int.sim.low <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.low[, i] <- tanh(dat_plot$B_cpc.1.3. + dat_plot$B_cpc.5.3. * (x2.sim[i]) + dat_plot$B_cpc.6.3. * (x2.sim[i]) * -1) 
}

# high density (+1SD)
int.sim.high <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim.high[, i] <- tanh(dat_plot$B_cpc.1.3. + dat_plot$B_cpc.5.3. * (x2.sim[i]) + dat_plot$B_cpc.6.3. * (x2.sim[i]) * 1) 
}


# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.3.) + median(dat_plot$B_cpc.5.3.) * (x2.sim) + median(dat_plot$B_cpc.6.3.) * (x2.sim) * -1)
bayes.c.eff.lower <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.low, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.low <- data.frame(x2.sim, bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

bayes.c.eff.mean <- tanh(median(dat_plot$B_cpc.1.3.) + median(dat_plot$B_cpc.5.3.) * (x2.sim) + median(dat_plot$B_cpc.6.3.) * (x2.sim) * 1)
bayes.c.eff.lower <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim.high, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat.high <- data.frame(x2.sim, bayes.c.eff.mean,
                            bayes.c.eff.lower, bayes.c.eff.upper,
                            bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                            bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                            bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                            bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

q3 <- ggplot(plot.dat.low, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#009E73", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#009E73", alpha = 0.1)+
  geom_line(color = "#009E73", size = 1.8, alpha=0.6)+
  
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "#E69F00", alpha = 0.1)+
  geom_ribbon(data=plot.dat.high, aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "#E69F00", alpha = 0.1)+
  geom_line(data=plot.dat.high, aes(x = x2.sim, y = bayes.c.eff.mean), color = "#E69F00", size = 1.8, alpha=0.6)+
  
  ylim(-1,1)+
  ggtitle("Brood size - Recruitment")+
  xlab("Beech mast index")+
  ylab("Observation-level correlation")+
  scale_x_continuous(breaks=c(0,1,2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
q3


#######################
#######################
#### Combine plots ####
#######################
#######################


combined_plot <-  (ann | ann_2) / (p | q) / (p2 | q2) / (p3 | q3) +
  plot_layout(heights = c(1, 8, 8, 8))
combined_plot

#ragg::agg_tiff("Figure S9.tiff", width = 8, height = 10, units = "in", res = 300)
combined_plot
#dev.off()

