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

load("~/data_wytham_tits.RData")
# this contains 4 objects:
# df: data frame containing data for the offspring mass model
# df_fecundity: data frame containing data for the fecundity and recruitment models
# mi_df: list of 20 data frames, with same format as "df" but with missing data imputed
# mi_df_fecundity: list of 20 data frames, with same format as "df_fecundity" but with missing data imputed

# see README file for more information on the dataframes

# remove data with NA, keeping only complete cases
df <- df[!is.na(df$Half_fall_date),]
df_fecundity <- df_fecundity[!is.na(df_fecundity$Half_fall_date),]
df <- df[!is.na(df$mast.score),]
df_fecundity <- df_fecundity[!is.na(df_fecundity$mast.score),]
df <- df[!is.na(df$breeding_age),]
df_fecundity <- df_fecundity[!is.na(df_fecundity$breeding_age),]
df <- df[!is.na(df$Mass.y),]
df_fecundity <- df_fecundity[!is.na(df_fecundity$Mass),]


##################
#### Analysis ####
##################
  
X <- as.matrix(cbind(rep(1, length(unique(df$BreedingSeason))),
                            scale(unique(df[,c("BreedingSeason","spring_temperature")])$spring_temperature)[,1],
                            scale(unique(df[,c("BreedingSeason","spring_precipitation")])$spring_precipitation)[,1],
                            scale(unique(df[,c("BreedingSeason","population_density")])$population_density)[,1],
                            unique(df[,c("BreedingSeason","mast.score")])$mast.score
  ))
  
  
X_g <- as.matrix(cbind(rep(1, dim(df)[1]),
                             scale(df$spring_temperature)[,1],
                             scale(df$spring_precipitation)[,1],
                             scale(df$population_density)[,1],
                             scale(df$Mass.y)[,1],
                             df$breeding_age,
                             scale(df$synchrony)[,1],
                             (scale(df$synchrony)[,1])^2,
                             df$mast.score
  ))
  
X_f <- as.matrix(cbind(
    scale(df_fecundity$spring_temperature)[,1],
    scale(df_fecundity$spring_precipitation)[,1],
    scale(df_fecundity$population_density)[,1],
    scale(df_fecundity$Mass)[,1],
    df_fecundity$breeding_age,
    scale(df_fecundity$synchrony)[,1],
    (scale(df_fecundity$synchrony)[,1])^2,
    df_fecundity$mast.score
  ))
  
X_r <- as.matrix(cbind(rep(1, dim(df_fecundity)[1]),
                             scale(df_fecundity$spring_temperature)[,1],
                             scale(df_fecundity$spring_precipitation)[,1],
                             scale(df_fecundity$population_density)[,1],
                             scale(df_fecundity$Mass)[,1],
                             df_fecundity$breeding_age,
                             scale(df_fecundity$synchrony)[,1],
                             (scale(df_fecundity$synchrony)[,1])^2,
                             df_fecundity$mast.score
  ))


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
rownames(df) <- NULL

  
stan.df =
    list(N = nrow(df),                      # total number of observations _ growth
         M = nrow(df_fecundity),            # total number of observations _ fecundity
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
         
         id_f = as.numeric(as.factor(df_fecundity$FemaleID)),  # index linking observations to individuals - fec
         c_id_f = as.numeric(as.factor(df_fecundity$BreedingSeason)), #index linking observations to contexts - fec
         idc_f = idc_f,                              # index linking individuals to positions in cmat - fecundity
         id_f_lm = as.numeric(rownames(df_fecundity)), # index observations - fecundity
         plot_id_f = as.numeric(as.factor(df_fecundity$LocationID)),
         
         X = X,                              # environmental predictor matrix (+ intercept) on correlation
         X_g = X_g,                          # environmental predictor matrix (+ intercept) on growth
         X_f = X_f,                          # environmental predictor matrix (+ intercept) on fecundity
         X_r = X_r,                          # environmental predictor matrix (+ intercept) on recruitment
         A = diag(length(unique(df$FemaleID))),   # relatedness matrix (identity matrix in this case)
         
         cm = max(as.numeric(table(as.numeric(as.factor(df_fecundity$BreedingSeason))))),
         cmat = cmat,
         cn = cn,
         cnt = length(unique(paste(df$FemaleID, df$BreedingSeason))), 
         
         
         N_id_zero = length(as.numeric(rownames(df_fecundity[df_fecundity$n_recruits==0,]))),
         N_id_nonzero = length(as.numeric(rownames(df_fecundity[df_fecundity$n_recruits!=0,]))),
         id_zero = as.numeric(rownames(df_fecundity[df_fecundity$n_recruits==0,])),
         id_nonzero = as.numeric(rownames(df_fecundity[df_fecundity$n_recruits!=0,])),
        
         
         growth = scale(as.numeric(df$Mass.x))[,1],   # offspring mass (response variable)
         productivity = df_fecundity$BroodSize_observed,    # fecundity (response variable))
         recruitment = df_fecundity$n_recruits
    )


########################
#### Run stan model ####
########################

library(shinystan)
library(cmdstanr)

mod <- cmdstan_model("model_wytham_tits.stan"
                     , stanc_options = list("O1")
)


fit <- mod$sample(
    data = stan.df, 
    output_dir = "/data/great_tit",
    seed = 4545, 
    chains = 3, 
    parallel_chains = 3,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.95,
    refresh = 20 # print update every 20 iters
  )

abcd <- fit$summary(variables = c("B_m_g", "B_m_f", "B_m_r", "B_cpc", "B_v", "cutpoint",
                                       "sd_g", "sd_f", "sd_r", "sd_E", "theta", "Sigma_plot", "sd_nestbox"))

# Combine the output of all the models
CRN_combined_fit <- read_cmdstan_csv(c(fit$output_files()),
                                     variables = c("B_m_g", "B_m_f", "B_m_r", "B_cpc", "B_v", "cutpoint",
                                                   "sd_g", "sd_f", "sd_r", "sd_E", "theta", "Sigma_plot", "sd_nestbox"))

dat_plot <- posterior::as_draws_df(CRN_combined_fit$post_warmup_draws)  
write.table(dat_plot, "posterior_draws_no_imputation.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
mcmc_trace(dat_plot)


#####################################
#### posterior predictive checks ####
#####################################

# Clutch size

draws_f <- fit$draws(
  variables = "y_rep_f",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_f <- ppc_dens_overlay(
  stan.df$productivity,
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

draws_g <- fit$draws(
  variables = "y_rep_g",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_g <- ppc_dens_overlay(
  stan.df$growth,
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

draws_r <- fit$draws(
  variables = "y_rep_r",
  inc_warmup = FALSE,
  format = "matrix"
)

plot_ppc_r <- ppc_dens_overlay(
  stan.df$recruitment,
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


