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


###################
#### Load data ####
###################

load("~/data_wytham_tits.RData")
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
  
  # Covariates on the correlations and on traits variances
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
                             mi_df[[i]]$mast.score
  ))
  
  # Covariates on brood size
  X_f[[i]] = as.matrix(cbind(
    scale(df_fecundity$spring_temperature)[,1],
    scale(df_fecundity$spring_precipitation)[,1],
    scale(df_fecundity$population_density)[,1],
    scale(mi_df_fecundity[[i]]$Mass)[,1],
    mi_df_fecundity[[i]]$breeding_age,
    scale(mi_df_fecundity[[i]]$synchrony)[,1],
    (scale(mi_df_fecundity[[i]]$synchrony)[,1])^2,
    mi_df_fecundity[[i]]$mast.score
  ))
  
  # Covariates on recruitment
  X_r[[i]] = as.matrix(cbind(rep(1, dim(df_fecundity)[1]),
                             scale(df_fecundity$spring_temperature)[,1],
                             scale(df_fecundity$spring_precipitation)[,1],
                             scale(df_fecundity$population_density)[,1],
                             scale(mi_df_fecundity[[i]]$Mass)[,1],
                             mi_df_fecundity[[i]]$breeding_age,
                             scale(mi_df_fecundity[[i]]$synchrony)[,1],
                             (scale(mi_df_fecundity[[i]]$synchrony)[,1])^2,
                             mi_df_fecundity[[i]]$mast.score
  ))
  
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
         M = nrow(df_fecundity),            # total number of observations _ fecundity/recruitment
         C = length(unique(df$BreedingSeason)),     # total number of environmental contexts (years)
         P = length(unique(df$LocationID)),     # total number of plots
         I = length(unique(df$FemaleID)),        # total number of subjects
         D = 3,                             # total number of traits/dimensions
         P_y = 5,
         P_g = 9,
         P_f = 8,
         P_r = 9,
         
         id_g = as.numeric(as.factor(df$FemaleID)),          # index linking observations to individuals - growth
         c_id_g = as.numeric(as.factor(df$BreedingSeason)),     # index linking observations to contexts - growth
         idc_g = idc_g,                                 # index linking individuals to positions in cmat - growth
         id_g_lm = as.numeric(rownames(df)),            # index observation - growth
         plot_id_g = as.numeric(as.factor(df$LocationID)),
         
         id_f = as.numeric(as.factor(df_fecundity$FemaleID)),  # index linking observations to individuals - fecundity/recruitment
         c_id_f = as.numeric(as.factor(df_fecundity$BreedingSeason)), #index linking observations to contexts - fecundity/recruitment
         idc_f = idc_f,                              # index linking individuals to positions in cmat - fecundity/recruitment
         id_f_lm = as.numeric(rownames(df_fecundity)), # index observations - fecundity/recruitment
         plot_id_f = as.numeric(as.factor(df_fecundity$LocationID)),
         
         X = X[[i]],                              # environmental predictor matrix (+ intercept) on correlation and variances
         X_g = X_g[[i]],                          # environmental predictor matrix (+ intercept) on growth
         X_f = X_f[[i]],                          # environmental predictor matrix (no intercept) on fecundity
         X_r = X_r[[i]],                          # environmental predictor matrix (+ intercept) on recruitment
         A = diag(length(unique(df$FemaleID))),   # relatedness matrix (identity matrix in this case)
         
         cm = max(as.numeric(table(as.numeric(as.factor(df_fecundity$BreedingSeason))))),
         cmat = cmat,
         cn = cn,
         cnt = length(unique(paste(df$FemaleID, df$BreedingSeason))), 
         
         id_zero = as.numeric(rownames(df_fecundity[df_fecundity$n_recruits==0,])),
         id_nonzero = as.numeric(rownames(df_fecundity[df_fecundity$n_recruits!=0,])),
         
         growth = as.numeric(df$Mass.x),             # offspring mass (response variable)
         productivity = df_fecundity$BroodSize_observed,    # fecundity (response variable))
         recruitment = df_fecundity$n_recruits           # offspring recruitment (response variable))
    )
  
}


########################
#### Run stan model ####
########################

library(shinystan)
library(cmdstanr)

mod <- cmdstan_model("model_wytham_tits.stan"
                     , stanc_options = list("O1")
)


# Run model on 20 alternative datasets
 
fit <- list()

for(i in 1:20) {
  
  fit[[i]] <- mod$sample(
    data = stan.df[[i]], 
    output_dir = "/data/great_tit",
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

# check ouput summary, e.g. here for the first model
abcd <- fit[[1]]$summary(variables = c("B_m_g", "B_m_f", "B_m_r", "B_cpc", "B_v", "cutpoint",
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
                                     variables = c("B_m_g", "B_m_f", "B_m_r", "B_cpc", "B_v", "cutpoint",
                                                   "sd_g", "sd_f", "sd_r", "sd_E", "theta", "Sigma_plot", "sd_nestbox"))

dat_plot <- posterior::as_draws_df(CRN_combined_fit$post_warmup_draws)  

# save posterior samples
write.table(dat_plot, "posterior_draws.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#check traceplots
mcmc_trace(dat_plot)



#####################################
#### posterior predictive checks ####
#####################################

# Clutch size

draws_f <- fit[[1]]$draws(
  variables = "y_rep_f",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_f <- ppc_dens_overlay(
  stan.df[[1]]$productivity,
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

draws_g <- fit[[1]]$draws(
  variables = "y_rep_g",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_g <- ppc_dens_overlay(
  stan.df[[1]]$growth,
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

draws_r <- fit[[1]]$draws(
  variables = "y_rep_r",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_r <- ppc_dens_overlay(
  stan.df[[1]]$recruitment,
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

