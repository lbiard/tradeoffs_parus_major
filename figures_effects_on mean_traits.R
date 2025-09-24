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

dat_plot <- read.csv("~/posterior_draws.txt", row.names=1, sep="")

df <- read.delim("~/df.txt")
df_fecundity <- read.delim("~/df_fecundity.txt")


################################
################################
#### Effects on mean traits ####
################################
################################


###########################################################
#### Functions to process output of ordinal regression ####
###########################################################

ordinal_fun <- function(i, x2.sim, dat_plot, linear_part) {
  
  # calculate probability of each brood size
  p0 <- ((exp(dat_plot$cutpoint.1. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.1. - linear_part * (x2.sim[i]))))
  p1 <- ((exp(dat_plot$cutpoint.2. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.2. - linear_part * (x2.sim[i])))) - p0
  p2 <- ((exp(dat_plot$cutpoint.3. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.3. - linear_part * (x2.sim[i])))) - p0 - p1
  p3 <- ((exp(dat_plot$cutpoint.4. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.4. - linear_part * (x2.sim[i])))) - p0 - p1 - p2
  p4 <- ((exp(dat_plot$cutpoint.5. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.5. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3
  p5 <- ((exp(dat_plot$cutpoint.6. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.6. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4
  p6 <- ((exp(dat_plot$cutpoint.7. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.7. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5
  p7 <- ((exp(dat_plot$cutpoint.8. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.8. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6
  p8 <- ((exp(dat_plot$cutpoint.9. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.9. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7
  p9 <- ((exp(dat_plot$cutpoint.10. - linear_part * (x2.sim[i])))/
           (1+exp(dat_plot$cutpoint.10. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8
  p10 <- ((exp(dat_plot$cutpoint.11. - linear_part * (x2.sim[i])))/
            (1+exp(dat_plot$cutpoint.11. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9
  p11 <- ((exp(dat_plot$cutpoint.12. - linear_part * (x2.sim[i])))/
            (1+exp(dat_plot$cutpoint.12. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9 - p10
  p12 <- ((exp(dat_plot$cutpoint.13. - linear_part * (x2.sim[i])))/
            (1+exp(dat_plot$cutpoint.13. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9 - p10 - p11
  p13 <- ((exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i])))/
            (1+exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i])))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9 - p10 - p11 - p12
  p14 <- 1 - ((exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i])))/
                (1+exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i]))))
  
  #transform back into an estimate on the scale of brood size
  estimated_brood <- 0*p0 + 1*p1 + 2*p2 + 3*p3 + 4*p4 + 5*p5 +
    6*p6 + 7*p7 + 8*p8 + 9*p9 + 10*p10 + 11*p11 +
    12*p12 + 13*p13 + 14*p14
  
  return(estimated_brood)
  
}


ordinal_fun_quad <- function(i, x2.sim, dat_plot, linear_part, quad_part) {
  
  # calculate probability of each brood size
  p0 <- ((exp(dat_plot$cutpoint.1. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.1. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2))))
  p1 <- ((exp(dat_plot$cutpoint.2. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.2. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0
  p2 <- ((exp(dat_plot$cutpoint.3. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.3. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1
  p3 <- ((exp(dat_plot$cutpoint.4. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.4. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2
  p4 <- ((exp(dat_plot$cutpoint.5. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.5. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3
  p5 <- ((exp(dat_plot$cutpoint.6. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.6. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4
  p6 <- ((exp(dat_plot$cutpoint.7. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.7. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5
  p7 <- ((exp(dat_plot$cutpoint.8. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.8. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6
  p8 <- ((exp(dat_plot$cutpoint.9. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.9. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7
  p9 <- ((exp(dat_plot$cutpoint.10. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
           (1+exp(dat_plot$cutpoint.10. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8
  p10 <- ((exp(dat_plot$cutpoint.11. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
            (1+exp(dat_plot$cutpoint.11. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9
  p11 <- ((exp(dat_plot$cutpoint.12. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
            (1+exp(dat_plot$cutpoint.12. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9 - p10
  p12 <- ((exp(dat_plot$cutpoint.13. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
            (1+exp(dat_plot$cutpoint.13. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9 - p10 - p11
  p13 <- ((exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
            (1+exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))) - p0 - p1 - p2 - p3 - p4 - p5 - p6 - p7 - p8 - p9 - p10 - p11 - p12
  p14 <- 1 - ((exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2)))/
                (1+exp(dat_plot$cutpoint.14. - linear_part * (x2.sim[i]) - quad_part * ((x2.sim[i])^2))))
  
  #transform back into an estimate on the scale of brood size
  estimated_brood <- 0*p0 + 1*p1 + 2*p2 + 3*p3 + 4*p4 + 5*p5 +
    6*p6 + 7*p7 + 8*p8 + 9*p9 + 10*p10 + 11*p11 +
    12*p12 + 13*p13 + 14*p14
  
  return(estimated_brood)
  
}



############################
#### Spring temperature ####
############################

## offspring mass

x2.sim = seq(min(scale(df$spring_temperature)[,1]),
             max(scale(df$spring_temperature)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$B_m_g.1. + dat_plot$B_m_g.2. * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.2.) * (x2.sim) 
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
inv_fun_p <- function(x){(x*sd(df$spring_temperature))+mean(df$spring_temperature)}


p_g_temperature <- ggplot(plot.dat, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df, aes(x = spring_temperature, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 0.03, height = 0)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring temperature")+
  ylab("Offpsring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## brood size

x2.sim = seq(min(scale(df_fecundity$spring_temperature)[,1]),
             max(scale(df_fecundity$spring_temperature)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ordinal_fun(i, x2.sim, dat_plot, linear_part = dat_plot$B_m_f.1.)
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, median)
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
inv_fun_p2 <- function(x){(x*sd(df_fecundity$spring_temperature))+mean(df_fecundity$spring_temperature)}


p_f_temperature <- ggplot(plot.dat, aes(x = inv_fun_p2(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = spring_temperature, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 0.03, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring temperature")+
  ylab("Brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## recruitment

x2.sim = seq(min(scale(df_fecundity$spring_temperature)[,1]),
             max(scale(df_fecundity$spring_temperature)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.2. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.2.) * (x2.sim))
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
inv_fun_p3 <- function(x){(x*sd(df_fecundity$spring_temperature))+mean(df_fecundity$spring_temperature)}


p_r_temperature <- ggplot(plot.dat, aes(x = inv_fun_p3(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = spring_temperature, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 0.03, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring temperature")+
  ylab("Recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


##############################
#### Spring precipitation ####
##############################

## offspring mass

x2.sim = seq(min(scale(df$spring_precipitation)[,1]),
             max(scale(df$spring_precipitation)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$B_m_g.1. + dat_plot$B_m_g.3. * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.3.) * (x2.sim) 
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
inv_fun_precipitation <- function(x){(x*sd(df$spring_precipitation))+mean(df$spring_precipitation)}


p_g_precipitation <- ggplot(plot.dat, aes(x = inv_fun_precipitation(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df, aes(x = spring_precipitation, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 1, height = 0)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring precipitation")+
  ylab("Offpsring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## brood size

x2.sim = seq(min(scale(df_fecundity$spring_precipitation)[,1]),
             max(scale(df_fecundity$spring_precipitation)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ordinal_fun(i, x2.sim, dat_plot, linear_part = dat_plot$B_m_f.2.)
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, median)
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
inv_fun_precipitation2 <- function(x){(x*sd(df_fecundity$spring_precipitation))+mean(df_fecundity$spring_precipitation)}


p_f_precipitation <- ggplot(plot.dat, aes(x = inv_fun_precipitation2(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = spring_precipitation, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 1, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring precipitation")+
  ylab("Brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## recruitment

x2.sim = seq(min(scale(df_fecundity$spring_precipitation)[,1]),
             max(scale(df_fecundity$spring_precipitation)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.3. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.3.) * (x2.sim))
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
inv_fun_precipitation3 <- function(x){(x*sd(df_fecundity$spring_precipitation))+mean(df_fecundity$spring_precipitation)}


p_r_precipitation <- ggplot(plot.dat, aes(x = inv_fun_precipitation3(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = spring_precipitation, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 1, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Spring precipitation")+
  ylab("Recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


############################
#### Population density ####
############################

## offspring mass

x2.sim = seq(min(scale(df$population_density)[,1]),
             max(scale(df$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$B_m_g.1. + dat_plot$B_m_g.4. * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.4.) * (x2.sim) 
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
inv_fun_density <- function(x){(x*sd(df$population_density))+mean(df$population_density)}


p_g_density <- ggplot(plot.dat, aes(x = inv_fun_density(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df, aes(x = population_density, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 2, height = 0)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Population density")+
  ylab("Offpsring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## brood size

x2.sim = seq(min(scale(df_fecundity$population_density)[,1]),
             max(scale(df_fecundity$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ordinal_fun(i, x2.sim, dat_plot, linear_part = dat_plot$B_m_f.3.)
}

# calculate quantiles of predictions
bayes.c.eff.mean <-  apply(int.sim, 2, median)
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
inv_fun_density2 <- function(x){(x*sd(df_fecundity$population_density))+mean(df_fecundity$population_density)}


p_f_density <- ggplot(plot.dat, aes(x = inv_fun_density2(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = population_density, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 2, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Population density")+
  ylab("Brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## recruitment

x2.sim = seq(min(scale(df_fecundity$population_density)[,1]),
             max(scale(df_fecundity$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.4. * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.4.) * (x2.sim))
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
inv_fun_density3 <- function(x){(x*sd(df_fecundity$population_density))+mean(df_fecundity$population_density)}


p_r_density <- ggplot(plot.dat, aes(x = inv_fun_density3(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = population_density, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 2, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Population density")+
  ylab("Recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



###################
#### Synchrony ####
###################

## offspring mass

x2.sim = seq(min(scale(df$synchrony)[,1], na.rm = T),
             max(scale(df$synchrony)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$B_m_g.1. + dat_plot$B_m_g.7. * (x2.sim[i]) + dat_plot$B_m_g.8. * (x2.sim[i] ^ 2) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.7.) * (x2.sim) + median(dat_plot$B_m_g.8.) * (x2.sim ^ 2) 
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
inv_fun_synchrony <- function(x){(x*sd(df$synchrony, na.rm = T))+mean(df$synchrony, na.rm = T)}


p_g_synchrony <- ggplot(plot.dat, aes(x = inv_fun_synchrony(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df, aes(x = synchrony, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 0.2, height = 0)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Synchrony")+
  ylab("Offpsring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## brood size

x2.sim = seq(min(scale(df_fecundity$synchrony)[,1], na.rm = T),
             max(scale(df_fecundity$synchrony)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ordinal_fun_quad(i, x2.sim, dat_plot, linear_part = dat_plot$B_m_f.6., quad_part = dat_plot$B_m_f.7.)
}

# calculate quantiles of predictions
bayes.c.eff.mean <-  apply(int.sim, 2, median)
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
inv_fun_synchrony2 <- function(x){(x*sd(df_fecundity$synchrony, na.rm = T))+mean(df_fecundity$synchrony, na.rm = T)}


p_f_synchrony <- ggplot(plot.dat, aes(x = inv_fun_synchrony2(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = synchrony, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 0.2, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Synchrony")+
  ylab("Brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## recruitment

x2.sim = seq(min(scale(df_fecundity$synchrony)[,1], na.rm = T),
             max(scale(df_fecundity$synchrony)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.7. * (x2.sim[i]) + dat_plot$B_m_r.8. * (x2.sim[i] ^ 2))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.7.) * (x2.sim) + median(dat_plot$B_m_r.8.) * (x2.sim^2))
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
inv_fun_synchrony3 <- function(x){(x*sd(df_fecundity$synchrony, na.rm = T))+mean(df_fecundity$synchrony, na.rm = T)}


p_r_synchrony <- ggplot(plot.dat, aes(x = inv_fun_synchrony3(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = synchrony, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 0.2, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Synchrony")+
  ylab("Recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


##########################
#### Beach mast index ####
##########################

## offspring mass

x2.sim = seq(0,
             2,
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$B_m_g.1. + dat_plot$B_m_g.9. * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.9.) * (x2.sim) 
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


p_g_beechmast <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_jitter(data = df, aes(x = mast.score, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 0.2, height = 0)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Beech mast index")+
  ylab("Offpsring mass")+
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## brood size

x2.sim = seq(0,
             2,
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ordinal_fun(i, x2.sim, dat_plot, linear_part = dat_plot$B_m_f.8.)
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, median)
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


p_f_beechmast <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = mast.score, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 0.2, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Beech mast index")+
  ylab("Brood size")+
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## recruitment

x2.sim = seq(0,
             2,
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.9. * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.9.) * (x2.sim))
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


p_r_beechmast <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = mast.score, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 0.2, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Beech mast index")+
  ylab("Recruitment")+
  scale_x_continuous(breaks=c(0,1,2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))



######################
#### Breeding age ####
######################

df$breeding_age_plot <- df$breeding_age+1
df_fecundity$breeding_age_plot <- df_fecundity$breeding_age+1

## offspring mass

int.sim <-  dat_plot$B_m_g.1.
int.sim2 <-  dat_plot$B_m_g.1. + dat_plot$B_m_g.6. 


# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) 
bayes.c.eff.lower <- quantile(int.sim, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim, probs = c(0.55))
plot.dat <- data.frame(bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.6.) 
bayes.c.eff.lower <- quantile(int.sim2, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim2, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim2, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim2, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim2, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim2, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim2, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim2, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim2, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim2, probs = c(0.55))
plot.dat.new <- data.frame(bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

plot.dat <- rbind(plot.dat, plot.dat.new)
sp <- as.data.frame(matrix(c("Young", "Old"), nrow = 2, ncol = 1))
plot.dat <- cbind(plot.dat, sp)


p_g_breedingage <- ggplot(df, aes(x = breeding_age, y = Mass.x)) +
  geom_jitter(data = df, aes(x = breeding_age_plot, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 0.3, height = 0)+
  xlab("Breeding age")+
  ylab("Offpsring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[1], ymax = plot.dat$bayes.c.eff.upper.bis4[1], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[2], ymax = plot.dat$bayes.c.eff.upper.bis4[2], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[1], ymax = plot.dat$bayes.c.eff.upper.bis3[1], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[2], ymax = plot.dat$bayes.c.eff.upper.bis3[2], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[1], ymax = plot.dat$bayes.c.eff.upper.bis2[1], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[2], ymax = plot.dat$bayes.c.eff.upper.bis2[2], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis[1], ymax = plot.dat$bayes.c.eff.upper.bis[1], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis[2], ymax = plot.dat$bayes.c.eff.upper.bis[2], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower[1], ymax = plot.dat$bayes.c.eff.upper[1], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower[2], ymax = plot.dat$bayes.c.eff.upper[2], alpha=0.1, fill="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.mean[1], yend = plot.dat$bayes.c.eff.mean[1], size=1.8, alpha=0.6, colour="darkseagreen4")
p_g_breedingage <- p_g_breedingage + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.mean[2], yend = plot.dat$bayes.c.eff.mean[2], size=1.8, alpha=0.6, colour="darkseagreen4")
p_g_breedingage <- p_g_breedingage + scale_x_continuous(breaks = 1:2,
                                                        labels = c("Young", "Old"))

## brood size

int.sim <- ordinal_fun(1, 0, dat_plot, linear_part = dat_plot$B_m_f.5.)
int.sim2 <- ordinal_fun(1, 1, dat_plot, linear_part = dat_plot$B_m_f.5.)

# calculate quantiles of predictions
bayes.c.eff.mean <- median(int.sim)
bayes.c.eff.lower <- quantile(int.sim, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim, probs = c(0.55))
plot.dat <- data.frame(bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

bayes.c.eff.mean <- median(int.sim2)
bayes.c.eff.lower <- quantile(int.sim2, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim2, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim2, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim2, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim2, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim2, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim2, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim2, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim2, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim2, probs = c(0.55))
plot.dat.new <- data.frame(bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

plot.dat <- rbind(plot.dat, plot.dat.new)
sp <- as.data.frame(matrix(c("Young", "Old"), nrow = 2, ncol = 1))
plot.dat <- cbind(plot.dat, sp)


p_f_breedingage <- ggplot(df_fecundity, aes(x = breeding_age, y = BroodSize_observed)) +
  geom_jitter(data = df_fecundity, aes(x = breeding_age_plot, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 0.3, height = 0.2)+
  xlab("Breeding age")+
  ylab("Brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[1], ymax = plot.dat$bayes.c.eff.upper.bis4[1], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[2], ymax = plot.dat$bayes.c.eff.upper.bis4[2], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[1], ymax = plot.dat$bayes.c.eff.upper.bis3[1], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[2], ymax = plot.dat$bayes.c.eff.upper.bis3[2], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[1], ymax = plot.dat$bayes.c.eff.upper.bis2[1], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[2], ymax = plot.dat$bayes.c.eff.upper.bis2[2], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis[1], ymax = plot.dat$bayes.c.eff.upper.bis[1], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis[2], ymax = plot.dat$bayes.c.eff.upper.bis[2], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower[1], ymax = plot.dat$bayes.c.eff.upper[1], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower[2], ymax = plot.dat$bayes.c.eff.upper[2], alpha=0.1, fill="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.mean[1], yend = plot.dat$bayes.c.eff.mean[1], size=1.8, alpha=0.6, colour="darkseagreen4")
p_f_breedingage <- p_f_breedingage + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.mean[2], yend = plot.dat$bayes.c.eff.mean[2], size=1.8, alpha=0.6, colour="darkseagreen4")
p_f_breedingage <- p_f_breedingage + scale_x_continuous(breaks = 1:2,
                                                        labels = c("Young", "Old"))

## recruitment 

int.sim <-  exp(dat_plot$B_m_r.1.)
int.sim2 <-  exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.6.)


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.))
bayes.c.eff.lower <- quantile(int.sim, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim, probs = c(0.55))
plot.dat <- data.frame(bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.6.))
bayes.c.eff.lower <- quantile(int.sim2, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim2, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim2, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim2, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim2, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim2, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim2, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim2, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim2, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim2, probs = c(0.55))
plot.dat.new <- data.frame(bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

plot.dat <- rbind(plot.dat, plot.dat.new)
sp <- as.data.frame(matrix(c("Young", "Old"), nrow = 2, ncol = 1))
plot.dat <- cbind(plot.dat, sp)


p_r_breedingage <- ggplot(df_fecundity, aes(x = breeding_age, y = n.recruits)) +
  geom_jitter(data = df_fecundity, aes(x = breeding_age_plot, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 0.3, height = 0.2)+
  xlab("Breeding age")+
  ylab("Recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[1], ymax = plot.dat$bayes.c.eff.upper.bis4[1], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[2], ymax = plot.dat$bayes.c.eff.upper.bis4[2], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[1], ymax = plot.dat$bayes.c.eff.upper.bis3[1], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[2], ymax = plot.dat$bayes.c.eff.upper.bis3[2], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[1], ymax = plot.dat$bayes.c.eff.upper.bis2[1], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[2], ymax = plot.dat$bayes.c.eff.upper.bis2[2], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis[1], ymax = plot.dat$bayes.c.eff.upper.bis[1], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis[2], ymax = plot.dat$bayes.c.eff.upper.bis[2], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower[1], ymax = plot.dat$bayes.c.eff.upper[1], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower[2], ymax = plot.dat$bayes.c.eff.upper[2], alpha=0.1, fill="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.mean[1], yend = plot.dat$bayes.c.eff.mean[1], size=1.8, alpha=0.6, colour="darkseagreen4")
p_r_breedingage <- p_r_breedingage + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.mean[2], yend = plot.dat$bayes.c.eff.mean[2], size=1.8, alpha=0.6, colour="darkseagreen4")
p_r_breedingage <- p_r_breedingage + scale_x_continuous(breaks = 1:2,
                                                        labels = c("Young", "Old"))


#####################
#### Mother mass ####
#####################

## offspring mass

x2.sim = seq(min(scale(df$Mass.y)[,1], na.rm = T),
             max(scale(df$Mass.y)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$B_m_g.1. + dat_plot$B_m_g.5. * (x2.sim[i]) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$B_m_g.1.) + median(dat_plot$B_m_g.5.) * (x2.sim)
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
inv_fun_mass <- function(x){(x*sd(df$Mass.y, na.rm = T))+mean(df$Mass.y, na.rm = T)}


p_g_mass <- ggplot(plot.dat, aes(x = inv_fun_mass(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df, aes(x = Mass.y, y = Mass.x), inherit.aes = F, alpha = 0.01, width = 0, height = 0)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Parental mass")+
  ylab("Offpsring mass")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## brood size

x2.sim = seq(min(scale(df_fecundity$Mass)[,1], na.rm = T),
             max(scale(df_fecundity$Mass)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- ordinal_fun(i, x2.sim, dat_plot, linear_part = dat_plot$B_m_f.4.)
}

# calculate quantiles of predictions
bayes.c.eff.mean <- apply(int.sim, 2, median)
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
inv_fun_mass2 <- function(x){(x*sd(df_fecundity$Mass, na.rm = T))+mean(df_fecundity$Mass, na.rm = T)}


p_f_mass <- ggplot(plot.dat, aes(x = inv_fun_mass2(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = Mass, y = BroodSize_observed), inherit.aes = F, alpha = 0.03, width = 0, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Parental mass")+
  ylab("Brood size")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## recruitment

x2.sim = seq(min(scale(df_fecundity$Mass)[,1], na.rm = T),
             max(scale(df_fecundity$Mass)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$B_m_r.1. + dat_plot$B_m_r.5. * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$B_m_r.1.) + median(dat_plot$B_m_r.5.) * (x2.sim))
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
inv_fun_mass3 <- function(x){(x*sd(df_fecundity$Mass, na.rm = T))+mean(df_fecundity$Mass, na.rm = T)}


p_r_mass <- ggplot(plot.dat, aes(x = inv_fun_mass3(x2.sim), y = bayes.c.eff.mean)) +
  geom_jitter(data = df_fecundity, aes(x = Mass, y = n_recruits), inherit.aes = F, alpha = 0.03, width = 0, height = 0.2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  xlab("Parental mass")+
  ylab("Recruitment")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


########################
#### Assemble plots ####
########################


#ragg::agg_tiff("Figure S7.tiff", width = 17, height = 10, units = "in", res = 300)
(p_g_temperature | p_g_precipitation | p_g_density | p_g_beechmast | p_g_breedingage | p_g_mass | p_g_synchrony) /
  (p_f_temperature | p_f_precipitation | p_f_density | p_f_beechmast | p_f_breedingage | p_f_mass | p_f_synchrony) /
  (p_r_temperature | p_r_precipitation | p_r_density | p_r_beechmast | p_r_breedingage | p_r_mass | p_r_synchrony)
#dev.off()





#####################################
#### Forest plot primary effects ####
#####################################


# change the "30000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 30 = 30000)

niter=60000 

df.posteriors <- data_frame(Submodel = c(rep("Offspring mass", niter*8), rep("Brood size", niter*8), rep("Recruitment", niter*8))
                            , parameter = c(rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Synchrony", niter), rep("Synchrony^2", niter), rep("Breeding age", niter), rep("Parental mass", niter)
                                            , rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Synchrony", niter), rep("Synchrony^2", niter), rep("Breeding age", niter), rep("Parental mass", niter)
                                            , rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Synchrony", niter), rep("Synchrony^2", niter), rep("Breeding age", niter), rep("Parental mass", niter))
                            , Posterior = c(dat_plot$B_m_g.4., dat_plot$B_m_g.2., dat_plot$B_m_g.3., dat_plot$B_m_g.9., dat_plot$B_m_g.7., dat_plot$B_m_g.8., dat_plot$B_m_g.6., dat_plot$B_m_g.5.
                                            , dat_plot$B_m_f.3., dat_plot$B_m_f.1., dat_plot$B_m_f.2., dat_plot$B_m_f.8., dat_plot$B_m_f.6., dat_plot$B_m_f.7., dat_plot$B_m_f.5., dat_plot$B_m_f.4.
                                            , dat_plot$B_m_r.4., dat_plot$B_m_r.2., dat_plot$B_m_r.3., dat_plot$B_m_r.9., dat_plot$B_m_r.7., dat_plot$B_m_r.8., dat_plot$B_m_r.6., dat_plot$B_m_r.5.))



df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("Offspring mass", "Brood size", "Recruitment"))
df.posteriors$parameter <- factor(df.posteriors$parameter, levels = c("Beech mast index", "Spring precipitation", "Spring temperature", "Population density", "Synchrony^2", "Synchrony", "Breeding age", "Parental mass"))

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
plot_forest <-ggplot()+
  stat_halfeye(data = df.posteriors, aes(x = Posterior, y = parameter, fill = Submodel), color = NA, alpha = 0.2, position = position_dodge(width = dodge.width), normalize="xy", scale=0.55)+
  geom_point(data = df.summary, aes(x = mean, y = parameter, color = Submodel), size = 2, position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI89_lower, xmax = BCI89_upper, y = parameter, color = Submodel), size=0.6, linetype="solid", position = position_dodge(width = dodge.width))+
  geom_linerange(data = df.summary, aes(xmin = BCI50_lower, xmax = BCI50_upper, y = parameter, color = Submodel), size=1.5, linetype="solid", position = position_dodge(width = dodge.width))+
  scale_color_manual(values = c(colT1, colT2, colCov))+
  scale_fill_manual(values = c(colT1, colT2, colCov))+
  scale_alpha_manual(values = c(0.3,1))+
  geom_vline(xintercept = 0, linetype = "dashed", color = colPlot, size = 0.4) +
  scale_x_continuous(breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5))+
  coord_cartesian(xlim=c(-0.7,0.7), clip = "off")+
  ylab("")+
  theme_minimal()+
  theme(plot.margin = margin(20,5.5,5.5,5.5))+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )

#ragg::agg_tiff("Figure S8.tiff", width = 8, height = 7, units = "in", res = 300)
plot_forest
#dev.off()





