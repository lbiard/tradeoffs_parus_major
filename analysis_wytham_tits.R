###################################
#### clean working environment ####
###################################

rm(list = ls())

#######################
#### load packages ####
#######################

library(dplyr)
library(mice)
library(patchwork)
library(ggplot2)
library(bayesplot)
library(purrr)
library(dplyr)
library(shinystan)
library(cmdstanr)
library(ggdist)


##############@####
#### Load data ####
###############@###

load("~/My Drive/phd/chapter 3/wytham_wood/analysis/last_version/cleaned data and scripts/data_wytham_tits.RData")
# this contains 4 objects:
# df: data frame containing data for the offspring mass model
# df_fecundity: data frame containing data for the fecundity and recruitment models
# mi_df: list of 20 data frames, with same format as "df" but with missing data imputed
# mi_df_fecundity: list of 20 data frames, with same format as "df_fecundity" but with missing data imputed

# see README file for more information on the dataframes


##################
#### Analysis ####
##################

X <- list()
X_g <- list()
X_f <- list()
X_r <- list()

for(i in 1:length(mi_df)) {

# Covariates on the correlations
X[[i]] <- as.matrix(cbind(rep(1, length(unique(df$BreedingSeason))),
                     scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1],
                     scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1],
                     scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1],
                     unique(mi_df[[i]][,c("BreedingSeason","mast.score")])$mast.score
                     ))

# Covariates on offspring mass
X_g[[i]] = as.matrix(cbind(rep(1, dim(df)[1]),
                      scale(df$spring_temperature)[,1],
                      scale(df$spring_precipitation)[,1],
                      scale(df$population_density)[,1],
                      scale(mi_df[[i]]$Mass)[,1],
                      mi_df[[i]]$breeding_age,
                      scale(mi_df[[i]]$synchrony)[,1],
                      (scale(mi_df[[i]]$synchrony)[,1])^2,
                      mi_df[[i]]$mast.score))

# Covariates on brood size
X_f[[i]] = as.matrix(cbind(rep(1, dim(df_fecundity)[1]),
                      scale(df_fecundity$spring_temperature)[,1],
                      scale(df_fecundity$spring_precipitation)[,1],
                      scale(df_fecundity$population_density)[,1],
                      scale(mi_df_fecundity[[i]]$Mass)[,1],
                      mi_df_fecundity[[i]]$breeding_age,
                      scale(mi_df_fecundity[[i]]$synchrony)[,1],
                      (scale(mi_df_fecundity[[i]]$synchrony)[,1])^2,
                      mi_df_fecundity[[i]]$mast.score))

# Covariates on recruitment
X_r[[i]] = as.matrix(cbind(rep(1, dim(df_fecundity)[1]),
                      scale(df_fecundity$spring_temperature)[,1],
                      scale(df_fecundity$spring_precipitation)[,1],
                      scale(df_fecundity$population_density)[,1],
                      scale(mi_df_fecundity[[i]]$Mass)[,1],
                      mi_df_fecundity[[i]]$breeding_age,
                      scale(mi_df_fecundity[[i]]$synchrony)[,1],
                      (scale(mi_df_fecundity[[i]]$synchrony)[,1])^2,
                      mi_df_fecundity[[i]]$mast.score))

}

cn <- as.numeric(table(as.numeric(as.factor(df_fecundity$BreedingSeason))))

#create empty cmat
cmat <- matrix(NA,
               nrow = length(unique(df$BreedingSeason)),
               ncol =  max(as.numeric(table(as.numeric(as.factor(df_fecundity$BreedingSeason))))))

#fill cmat
temporary <- as.data.frame(cbind(as.numeric(as.factor(df_fecundity$BreedingSeason)),
                                 as.numeric(as.factor(df_fecundity$FemaleID))))
for (i in 1:length(unique(df$BreedingSeason))) {
  cmat[i, 1:cn[i]] <- temporary$V2[temporary$V1 == i]
}
cmat_n = apply(cmat, 1, FUN = function(x) sum(!is.na(x)) )
cmat[is.na(cmat)] = 0 #remove NAs

temp = t(cmat)
corder = data.frame(id = temp[temp>0], c = rep(seq(1:nrow(cmat)), times = cmat_n))
idc_f = match(paste0(as.numeric(as.factor(df_fecundity$FemaleID)),
                     sep=".",
                     as.numeric(as.factor(df_fecundity$BreedingSeason))),
              paste0(corder$id, sep=".", corder$c))

idc_g = match(paste0(as.numeric(as.factor(df$FemaleID)),
                     sep=".",
                     as.numeric(as.factor(df$BreedingSeason))),
              paste0(corder$id, sep=".", corder$c))

rownames(df_fecundity) <- NULL

stan.df <- list()

for(i in 1:length(mi_df)) {
  
stan.df[[i]] =
  list(N = nrow(df),                      # total number of observations _ growth
       M = nrow(df_fecundity),            # total number of observations _ fecundity
       C = length(unique(df$BreedingSeason)),     # total number of environmental contexts (years)
       P = length(unique(df$LocationID)),     # total number of plots
       I = length(unique(df$FemaleID)),        # total number of subjects
       D = 3,                             # total number of traits/dimensions
       P_y = 5,       # number of covariates on correlations
       P_g = 9,       # number of covariates on offpsring mass
       P_f = 9,       # number of covariates on brood size
       P_r = 9,       # number of covariates on recruitment

       id_g = as.numeric(as.factor(df$FemaleID)),          # index linking observations to individuals - growth
       c_id_g = as.numeric(as.factor(df$BreedingSeason)),     # index linking observations to contexts - growth
       idc_g = idc_g,                                 # index linking individuals to positions in cmat - growth
       id_g_lm = as.numeric(rownames(df)),            # index observation - growth
       plot_id_g = as.numeric(as.factor(df$LocationID)),

       id_f = as.numeric(as.factor(df_fecundity$FemaleID)),  # index linking observations to individuals - fec
       c_id_f = as.numeric(as.factor(df_fecundity$BreedingSeason)), #index linking observations to contexts - fec
       idc_f = idc_f,                              # index linking individuals to positions in cmat - fecundity
       id_f_lm = as.numeric(rownames(df_fecundity)), # index observations - fecundity
       plot_id_f = as.numeric(as.factor(df_fecundity$LocationID)),

       X = X[[i]],                              # environmental predictor matrix (+ intercept) on correlation
       X_g = X_g[[i]],                          # environmental predictor matrix (+ intercept) on growth
       X_f = X_f[[i]],                          # environmental predictor matrix (+ intercept) on fecundity
       X_r = X_r[[i]],                          # environmental predictor matrix (+ intercept) on recruitment
       A = diag(length(unique(df$FemaleID))),   # relatedness matrix (identity matrix in this case)

       cm = max(as.numeric(table(as.numeric(as.factor(df_fecundity$BreedingSeason))))),
       cmat = cmat,
       cn = cn,
       cnt = length(unique(paste(df$FemaleID, df$BreedingSeason))), 

       id_zero = as.numeric(rownames(df_fecundity[df_fecundity$n_recruits==0,])),
       id_nonzero = as.numeric(rownames(df_fecundity[df_fecundity$n_recruits!=0,])),

       growth = as.numeric(df$Mass.x),                    # offspring mass (response variable)
       productivity = df_fecundity$BroodSize_observed,    # fecundity (response variable))
       recruitment = df_fecundity$n_recruits              # recruitment (response variable)
  )

}

########################
#### Run stan model ####
########################

setwd("~/My Drive/phd/chapter 3/wytham_wood/analysis/last_version/cleaned data and scripts")

library(shinystan)
library(cmdstanr)

mod <- cmdstan_model("model_wytham_tits.stan"
                     , stanc_options = list("O1")
)

fit <- list()

for(i in 1:20) {

fit[[i]] <- mod$sample(
  data = stan.df[[i]], 
  output_dir = "/data/lbliar/great_tit",
  seed = 125, 
  chains = 3, 
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  refresh = 20 # print update every 20 iters
)

print(i)

}

# check summary of the first model
abcd <- fit[[1]]$summary(variables = c("B_m_g", "B_m_f", "B_m_r", "B_cpc",
                               "sd_g", "sd_f", "sd_r", "sd_E", "theta", "Sigma_plot", "sd_nestbox"))

# Combine the output of all the models
CRN_combined_fit <- read_cmdstan_csv(c(fit[[1]]$output_files(), 
                                       fit[[2]]$output_files(),
                                       fit[[3]]$output_files(),
                                       fit[[4]]$output_files(),
                                       fit[[5]]$output_files(),
                                       fit[[6]]$output_files(), 
                                       fit[[7]]$output_files(),
                                       fit[[8]]$output_files(),
                                       fit[[9]]$output_files(),
                                       fit[[10]]$output_files(),
                                       fit[[11]]$output_files(), 
                                       fit[[12]]$output_files(),
                                       fit[[13]]$output_files(),
                                       fit[[14]]$output_files(),
                                       fit[[15]]$output_files(),
                                       fit[[16]]$output_files(), 
                                       fit[[17]]$output_files(),
                                       fit[[18]]$output_files(),
                                       fit[[19]]$output_files(),
                                       fit[[20]]$output_files()),
                         variables = c("B_m_g", "B_m_f", "B_m_r", "B_cpc",
                                       "sd_g", "sd_f", "sd_r", "sd_E", "theta", "Sigma_plot", "sd_nestbox"))

# Get pooled posterior draws
dat_plot <- posterior::as_draws_df(CRN_combined_fit$post_warmup_draws)  
# write.table(dat_plot, "posterior_draws.txt", append = FALSE, sep = " ", dec = ".",
#             row.names = TRUE, col.names = TRUE)

# check traceplots
mcmc_trace(dat_plot)


#####################################
#### posterior predictive checks ####
#####################################

# Clutch size

draws_f <- fit[[5]]$draws(
      variables = "y_rep_f",
      inc_warmup = FALSE,
      format = "matrix"
  )

plot_ppc_f <- ppc_dens_overlay(
      stan.df[[5]]$productivity,
      draws_f[sample(nrow(draws_f), 100),],
      size = 0.25,
      alpha = 0.7,
      trim = FALSE
  )

plot_ppc_f <- plot_ppc_f + 
      xlab("Clutch size") + 
      ylab("Density") + 
      ggtitle("Posterior predictive check \u2014 Clutch size")
plot_ppc_f

# Offspring mass

draws_g <- fit[[5]]$draws(
  variables = "y_rep_g",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_g <- ppc_dens_overlay(
  stan.df[[5]]$growth,
  draws_g[sample(nrow(draws_g), 100),],
  size = 0.25,
  alpha = 0.7,
  trim = FALSE
)

plot_ppc_g <- plot_ppc_g + 
  xlab("Offspring mass") + 
  ylab("Density") + 
  ggtitle("Posterior predictive check \u2014 Offspring mass")
plot_ppc_g

# Recruitment

draws_r <- fit[[5]]$draws(
  variables = "y_rep_r",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_r <- ppc_dens_overlay(
  stan.df[[5]]$recruitment,
  draws_r[sample(nrow(draws_r), 100),],
  size = 0.25,
  alpha = 0.7,
  trim = FALSE
)

plot_ppc_r <- plot_ppc_r + 
  xlab("Recruitment") + 
  ylab("Density") + 
  xlim(c(0,10))+
  ggtitle("Posterior predictive check \u2014 Recruitment")
plot_ppc_r

wrap_elements(full= (plot_ppc_f / plot_ppc_g / plot_ppc_r))

# remove large elements from memory
rm(draws_f, draws_g, draws_r)


#################
#################
#### Figures ####
#################
#################

############################
#### population density ####
############################

## Offpsring mass - Brood size ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
          max(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
          by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,1]` + dat_plot$`B_cpc[4,1]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,1]`) + median(dat_plot$`B_cpc[4,1]`) * (x2.sim)) 
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

inv_fun_p <- function(x){(x*sd(unique(df[,c("BreedingSeason","population_density")])$population_density))+mean(unique(df[,c("BreedingSeason","population_density")])$population_density)}


p <- ggplot(plot.dat, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Offpsring mass - Brood size")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

## Offpsring mass - Recruitment ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
          max(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
          by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,2]` + dat_plot$`B_cpc[4,2]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,2]`) + median(dat_plot$`B_cpc[4,2]`) * (x2.sim)) 
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

inv_fun_p2 <- function(x){(x*sd(unique(df[,c("BreedingSeason","population_density")])$population_density))+mean(unique(df[,c("BreedingSeason","population_density")])$population_density)}


q <- ggplot(plot.dat, aes(x = inv_fun_p2(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Offpsring mass - Recruitment")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## Brood size - Recruitment ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
          max(scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1]),
          by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,3]` + dat_plot$`B_cpc[4,3]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,3]`) + median(dat_plot$`B_cpc[4,3]`) * (x2.sim)) 
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

inv_fun_p3 <- function(x){(x*sd(unique(df[,c("BreedingSeason","population_density")])$population_density))+mean(unique(df[,c("BreedingSeason","population_density")])$population_density)}


r <- ggplot(plot.dat, aes(x = inv_fun_p3(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Brood size - Recruitment")+
  xlab("Population density")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

combined_plot <- p+q+r & xlab(NULL) & ylab(NULL)


#ragg::agg_tiff("Figure 1.tiff", width = 10, height = 4.5, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
abcd <- gridExtra::grid.arrange(gt, left = "Correlation", bottom = "Population density")  
#dev.off()



############################
#### Spring temperature ####
############################

## Offpsring mass - Brood size ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,1]` + dat_plot$`B_cpc[2,1]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,1]`) + median(dat_plot$`B_cpc[2,1]`) * (x2.sim)) 
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
inv_fun_p <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature))+mean(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)}


p <- ggplot(plot.dat, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Offpsring mass - Brood size")+
  xlab("Spring temperature")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

## Offpsring mass - Recruitment ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,2]` + dat_plot$`B_cpc[2,2]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,2]`) + median(dat_plot$`B_cpc[2,2]`) * (x2.sim)) 
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

inv_fun_p2 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature))+mean(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)}


q <- ggplot(plot.dat, aes(x = inv_fun_p2(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Offpsring mass - Recruitment")+
  xlab("Spring temperature")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## Brood size - Recruitment ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,3]` + dat_plot$`B_cpc[2,3]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,3]`) + median(dat_plot$`B_cpc[2,3]`) * (x2.sim)) 
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

inv_fun_p3 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature))+mean(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)}


r <- ggplot(plot.dat, aes(x = inv_fun_p3(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Brood size - Recruitment")+
  xlab("Spring temperature")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

combined_plot <- p+q+r & xlab(NULL) & ylab(NULL)

#ragg::agg_tiff("Figure 2.tiff", width = 10, height = 4.5, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
gridExtra::grid.arrange(gt, left = "Correlation", bottom = "Spring temperature")  
#dev.off()



##############################
#### Spring precipitation ####
##############################

## Offpsring mass - Brood size ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,1]` + dat_plot$`B_cpc[3,1]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,1]`) + median(dat_plot$`B_cpc[3,1]`) * (x2.sim)) 
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

inv_fun_p <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation))+mean(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)}


p <- ggplot(plot.dat, aes(x = inv_fun_p(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Offpsring mass - Brood size")+
  xlab("Spring precipitation")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

## Offpsring mass - Recruitment ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,2]` + dat_plot$`B_cpc[3,2]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,2]`) + median(dat_plot$`B_cpc[3,2]`) * (x2.sim)) 
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
inv_fun_p2 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation))+mean(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)}


q <- ggplot(plot.dat, aes(x = inv_fun_p2(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Offpsring mass - Recruitment")+
  xlab("Spring precipitation")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))


## Brood size - Recruitment ##

x2.sim = seq(min(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             max(scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1]),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,3]` + dat_plot$`B_cpc[3,3]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,3]`) + median(dat_plot$`B_cpc[3,3]`) * (x2.sim)) 
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

inv_fun_p3 <- function(x){(x*sd(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation))+mean(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)}


r <- ggplot(plot.dat, aes(x = inv_fun_p3(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
  ylim(-1,1)+
  ggtitle("Brood size - Recruitment")+
  xlab("Spring precipitation")+
  ylab("Observation-level correlation")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14))

combined_plot <- p+q+r & xlab(NULL) & ylab(NULL)

#ragg::agg_tiff("Figure 3.tiff", width = 10, height = 4.5, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
gridExtra::grid.arrange(gt, left = "Correlation", bottom = "Spring precipitation")  
#dev.off()


####################
#### Beech mast ####
####################

## Offpsring mass - Brood size ##

x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,1]` + dat_plot$`B_cpc[5,1]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,1]`) + median(dat_plot$`B_cpc[5,1]`) * (x2.sim)) 
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


p <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
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


## Offpsring mass - Recruitment ##

x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,2]` + dat_plot$`B_cpc[5,2]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,2]`) + median(dat_plot$`B_cpc[5,2]`) * (x2.sim)) 
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



q <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
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


## Brood size - Recruitment ##

x2.sim = seq(min(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             max(unique(df[,c("BreedingSeason","mast.score")])$mast.score, na.rm = T),
             by =  0.2) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- tanh(dat_plot$`B_cpc[1,3]` + dat_plot$`B_cpc[5,3]` * (x2.sim[i])) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- tanh(median(dat_plot$`B_cpc[1,3]`) + median(dat_plot$`B_cpc[5,3]`) * (x2.sim)) 
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


r <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "darkseagreen4", alpha = 0.1)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "darkseagreen4", alpha = 0.1)+
  geom_line(color = "darkseagreen4", size = 1.8, alpha=0.6)+
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

combined_plot <- p+q+r & xlab(NULL) & ylab(NULL)

#ragg::agg_tiff("Figure 4.tiff", width = 10, height = 4.5, units = "in", res = 300)
gt <- patchwork::patchworkGrob(combined_plot)
gridExtra::grid.arrange(gt, left = "Correlation", bottom = "Beech mast index")  
#dev.off()




###########################
#### Figure posteriors ####
###########################


# change the "60000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 60 = 60000)

niter=60000 

df.posteriors <- data_frame(Submodel = c(rep("offspring mass - brood size", niter*5), rep("offpsring mass - recruitment", niter*5), rep("brood size - recruitment", niter*5))
                            , parameter = c(rep("Intercept", niter), rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter)
                                            , rep("Intercept", niter), rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter)
                                            , rep("Intercept", niter), rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter))
                            , Posterior = c(dat_plot$`B_cpc[1,1]`, dat_plot$`B_cpc[4,1]`, dat_plot$`B_cpc[2,1]`, dat_plot$`B_cpc[3,1]`, dat_plot$`B_cpc[5,1]`
                                            , dat_plot$`B_cpc[1,2]`, dat_plot$`B_cpc[4,2]`, dat_plot$`B_cpc[2,2]`, dat_plot$`B_cpc[3,2]`, dat_plot$`B_cpc[5,2]`
                                            , dat_plot$`B_cpc[1,3]`, dat_plot$`B_cpc[4,3]`, dat_plot$`B_cpc[2,3]`, dat_plot$`B_cpc[3,3]`, dat_plot$`B_cpc[5,3]`))



df.posteriors$Submodel <- factor(df.posteriors$Submodel, levels=c("offspring mass - brood size", "offpsring mass - recruitment", "brood size - recruitment"))
df.posteriors$parameter <- factor(df.posteriors$parameter, levels = c("Beech mast index", "Spring precipitation", "Spring temperature", "Population density", "Intercept"))

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
p2 <-ggplot()+
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

#ragg::agg_tiff("Figure 5.tiff", width = 8, height = 7, units = "in", res = 300)
p2
#dev.off()


################################
#### Effects on mean traits ####
################################


#### Spring temperature ####

## offspring mass

x2.sim = seq(min(scale(df$spring_temperature)[,1]),
             max(scale(df$spring_temperature)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[2]` * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[2]`) * (x2.sim) 
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
  int.sim[, i] <- exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[2]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[2]`) * (x2.sim))
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
  int.sim[, i] <- exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[2]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[2]`) * (x2.sim))
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



#### Spring precipitation ####

## offspring mass

x2.sim = seq(min(scale(df$spring_precipitation)[,1]),
             max(scale(df$spring_precipitation)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[3]` * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[3]`) * (x2.sim) 
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
  int.sim[, i] <- exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[3]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[3]`) * (x2.sim))
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
  int.sim[, i] <- exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[3]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[3]`) * (x2.sim))
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



#### Population density ####

## offspring mass

x2.sim = seq(min(scale(df$population_density)[,1]),
             max(scale(df$population_density)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[4]` * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[4]`) * (x2.sim) 
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
  int.sim[, i] <- exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[4]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[4]`) * (x2.sim))
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
  int.sim[, i] <- exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[4]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[4]`) * (x2.sim))
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




#### Synchrony ####

## offspring mass

x2.sim = seq(min(scale(df$synchrony)[,1], na.rm = T),
             max(scale(df$synchrony)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[7]` * (x2.sim[i]) + dat_plot$`B_m_g[8]` * (x2.sim[i] ^ 2) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[7]`) * (x2.sim) + median(dat_plot$`B_m_g[8]`) * (x2.sim ^ 2) 
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
  int.sim[, i] <- exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[7]` * (x2.sim[i]) + dat_plot$`B_m_f[8]` * (x2.sim[i] ^ 2))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[7]`) * (x2.sim) + median(dat_plot$`B_m_f[8]`) * (x2.sim^2))
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
  int.sim[, i] <- exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[7]` * (x2.sim[i]) + dat_plot$`B_m_r[8]` * (x2.sim[i] ^ 2))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[7]`) * (x2.sim) + median(dat_plot$`B_m_r[8]`) * (x2.sim^2))
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



#### Beach mast index ####

## offspring mass

x2.sim = seq(0,
             2,
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[9]` * (x2.sim[i])
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[9]`) * (x2.sim) 
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
  int.sim[, i] <- exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[9]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[9]`) * (x2.sim))
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
  int.sim[, i] <- exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[9]` * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[9]`) * (x2.sim))
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




#### Breeding age ####

## offspring mass

int.sim <-  dat_plot$`B_m_g[1]`
int.sim2 <-  dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[6]` 


# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) 
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

bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[6]`) 
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

int.sim <-  exp(dat_plot$`B_m_f[1]`)
int.sim2 <-  exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[6]`)


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`))
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

bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[6]`))
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

int.sim <-  exp(dat_plot$`B_m_r[1]`)
int.sim2 <-  exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[6]`)


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`))
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

bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[6]`))
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



#### Mother mass ####

## offspring mass

x2.sim = seq(min(scale(df$Mass.y)[,1], na.rm = T),
             max(scale(df$Mass.y)[,1], na.rm = T),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- dat_plot$`B_m_g[1]` + dat_plot$`B_m_g[5]` * (x2.sim[i]) 
}

# calculate quantiles of predictions
bayes.c.eff.mean <- median(dat_plot$`B_m_g[1]`) + median(dat_plot$`B_m_g[5]`) * (x2.sim)
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
  int.sim[, i] <- exp(dat_plot$`B_m_f[1]` + dat_plot$`B_m_f[5]` * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_f[1]`) + median(dat_plot$`B_m_f[5]`) * (x2.sim))
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
  int.sim[, i] <- exp(dat_plot$`B_m_r[1]` + dat_plot$`B_m_r[5]` * (x2.sim[i]))
}


# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$`B_m_r[1]`) + median(dat_plot$`B_m_r[5]`) * (x2.sim))
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


#ragg::agg_tiff("Figure supplement.tiff", width = 14, height = 8, units = "in", res = 300)
(p_g_temperature | p_g_precipitation | p_g_density | p_g_beechmast | p_g_breedingage | p_g_mass | p_g_synchrony) /
  (p_f_temperature | p_f_precipitation | p_f_density | p_f_beechmast | p_f_breedingage | p_f_mass | p_f_synchrony) /
  (p_r_temperature | p_r_precipitation | p_r_density | p_r_beechmast | p_r_breedingage | p_r_mass | p_r_synchrony)
#dev.off()





#####################################
#### Forest plot primary effects ####
#####################################


# change the "60000" to another value if you run chains for less or more iterations than I did (1000 iterations post burn-in * 60 = 60000)

niter=60000 

df.posteriors <- data_frame(Submodel = c(rep("Offspring mass", niter*8), rep("Brood size", niter*8), rep("Recruitment", niter*8))
                            , parameter = c(rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Synchrony", niter), rep("Synchrony^2", niter), rep("Breeding age", niter), rep("Parental mass", niter)
                                            , rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Synchrony", niter), rep("Synchrony^2", niter), rep("Breeding age", niter), rep("Parental mass", niter)
                                            , rep("Population density", niter), rep("Spring temperature", niter), rep("Spring precipitation", niter), rep("Beech mast index", niter), rep("Synchrony", niter), rep("Synchrony^2", niter), rep("Breeding age", niter), rep("Parental mass", niter))
                            , Posterior = c(dat_plot$`B_m_g[4]`, dat_plot$`B_m_g[2]`, dat_plot$`B_m_g[3]`, dat_plot$`B_m_g[9]`, dat_plot$`B_m_g[7]`, dat_plot$`B_m_g[8]`, dat_plot$`B_m_g[6]`, dat_plot$`B_m_g[5]`
                                            , dat_plot$`B_m_f[4]`, dat_plot$`B_m_f[2]`, dat_plot$`B_m_f[3]`, dat_plot$`B_m_f[9]`, dat_plot$`B_m_f[7]`, dat_plot$`B_m_f[8]`, dat_plot$`B_m_f[6]`, dat_plot$`B_m_f[5]`
                                            , dat_plot$`B_m_r[4]`, dat_plot$`B_m_r[2]`, dat_plot$`B_m_r[3]`, dat_plot$`B_m_r[9]`, dat_plot$`B_m_r[7]`, dat_plot$`B_m_r[8]`, dat_plot$`B_m_r[6]`, dat_plot$`B_m_r[5]`))



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
  scale_x_continuous(breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4))+
  coord_cartesian(xlim=c(-0.5,0.5), clip = "off")+
  ylab("")+
  theme_minimal()+
  theme(plot.margin = margin(20,5.5,5.5,5.5))+
  theme(axis.line.y = element_blank()
        , axis.ticks.y = element_blank()
        , panel.grid.major = element_blank() 
        , panel.grid.minor = element_blank()
  )

#ragg::agg_tiff("Figure S2.tiff", width = 8, height = 7, units = "in", res = 300)
plot_forest
#dev.off()