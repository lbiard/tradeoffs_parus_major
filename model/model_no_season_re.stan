functions {

//functions are used fom prior work by
//Dan Schrage (https://gitlab.com/dschrage/rcovreg)
  
  real sum_square_x(matrix x, int i, int j) {
    int j_prime;
    real sum_x = 0;
    if(j==1) return(sum_x);

    j_prime = 1;
    while(j_prime < j) {
      sum_x += x[i,j_prime]^2;
      j_prime += 1;
    }
    return(sum_x);
  }
  
  matrix lkj_to_chol_corr(row_vector constrained_reals, int ntrait) {
    int z_counter;
    matrix[ntrait,ntrait] x;

    z_counter = 1;
    x[1,1] = 1;
    for(j in 2:ntrait) {
      x[1,j] = 0;
    }
    for(i in 2:ntrait) {
      for(j in 1:ntrait) {
        if(i==j) {
          x[i,j] = sqrt(1 - sum_square_x(x, i, j));
        } else if(i > j) {
          x[i,j] = constrained_reals[z_counter]*sqrt(1 - sum_square_x(x, i, j));
          z_counter += 1;
        } else {
          x[i,j] = 0;
        }
      }
    }
    return(x);
  }
  
  int num_zeros(array[] int y) {
    int sum = 0;
    for (m in 1:size(y)) {
      sum += (y[m] == 0);
    }
    return sum;
  }
  
}

data {
  int<lower=1> N; //total number of observations _ growth
  int<lower=1> M; //total number of observations - fecundity / recruitment
  int<lower=1> C; //total number of environmental contexts (years)
  int<lower=1> P; //total number of plots (subdivisions of Wytham)
  int<lower=1> I; //total number of subjects
  int<lower=1> D; //total number of traits/dimensions
  int<lower=0> P_y; //total number of environmental predictors (+ intercept) on correlation
  int<lower=0> P_g; //total number of environmental predictors (+ intercept) on growth
  int<lower=0> P_f; //total number of environmental predictors (+ intercept) on fecundity
  int<lower=0> P_r; //total number of environmental predictors (+ intercept) on recruitment
  
  array[N] int<lower=0> id_g; //index linking observations to individuals - growth
  array[N] int<lower=0> c_id_g; //index linking observations to contexts - growth
  array[N] int<lower=0> idc_g; //index linking individuals to positions in cmat - growth
  array[N] int<lower=0> id_g_lm; //index observation - growth
  array[N] int<lower=0> plot_id_g; //plot for each observation - growth
  
  array[M] int<lower=0> id_f; //index linking observations to individuals - fecundity
  array[M] int<lower=0> c_id_f; //index linking observations to contexts - fecundity
  array[M] int<lower=0> idc_f; //index linking individuals to positions in cmat - fecundity
  array[M] int<lower=0> id_f_lm; //index observations - fecundity
  array[M] int<lower=0> plot_id_f; //plot for each observation - fecundity
  
  matrix[C,P_y] X; //environmental predictor matrix (+ intercept) on correlation
  matrix[N,P_g] X_g; //environmental predictor matrix (+ intercept) on growth
  matrix[M,P_f] X_f; //environmental predictor matrix (+ intercept) on fecundity
  matrix[M,P_r] X_r; //environmental predictor matrix (+ intercept) on fecundity
  matrix[I,I] A; //relatedness matrix
  
  int<lower=1> cm; //max number of individuals observed in a context
  array[C, cm] int cmat; //matrix with all individuals observed in each context (row)
  array [C] int<lower=0> cn; //count of individuals observed per context
  int<lower=1> cnt; //total number of individuals across contexts
  
  int<lower=1> N_id_zero;   //number of zero recruitment
  int<lower=1> N_id_nonzero;  // number of non zero recruitment
  array[N_id_zero] int<lower=0> id_zero;
  array[N_id_nonzero] int<lower=0> id_nonzero;
  
  array[N] real growth; //offspring mass data
  array[M] int<lower=1, upper=15> productivity; //fecundity data
  array[M] int recruitment; //recruitment data
}

transformed data{
  matrix[I, I] LA = cholesky_decompose(A);
  int ncor = (D*(D-1))/2; //unique cov/cor parameters
  // Compute, thin, and then scale QR decomposition
  matrix[C, P_y] Q = qr_thin_Q(X) * sqrt(C-1);
  matrix[P_y, P_y] R = qr_thin_R(X) / sqrt(C-1);
  matrix[P_y, P_y] R_inv = inverse(R);
  
  matrix[N, P_g] Q_g = qr_thin_Q(X_g) * sqrt(N-1);
  matrix[P_g, P_g] R_g = qr_thin_R(X_g) / sqrt(N-1);
  matrix[P_g, P_g] R_inv_g = inverse(R_g);
  
  matrix[M, P_f] Q_f = qr_thin_Q(X_f) * sqrt(M-1);
  matrix[P_f, P_f] R_f = qr_thin_R(X_f) / sqrt(M-1);
  matrix[P_f, P_f] R_inv_f = inverse(R_f);
  
  matrix[M, P_r] Q_r = qr_thin_Q(X_r) * sqrt(M-1);
  matrix[P_r, P_r] R_r = qr_thin_R(X_r) / sqrt(M-1);
  matrix[P_r, P_r] R_inv_r = inverse(R_r);

  int<lower=0> M_zero = num_zeros(recruitment);
    array[M - M_zero] int<lower=1> y_nonzero;
    int M_nonzero = 0;
    for (m in 1:M) {
      if (recruitment[m] == 0) continue;
      M_nonzero += 1;
      y_nonzero[M_nonzero] = recruitment[m];
    }

}

parameters { 
  //fixed effects
  vector[P_g] B_mq_g; //RN of means - growth
  vector[P_f] B_mq_f; //RN of means - fecundity
  vector[P_r] B_mq_r; //RN of means - recruitment
  matrix[P_y, ncor] B_cpcq; //RN of canonical partial correlations
  matrix[P_y, D] B_vq; //RN of variances

  //random effects
  matrix[cnt, D] Z_G; //all context-specific additive genetic values
  real<lower=0> sd_E; //residual standard deviation (within litter variance) - growth
  // array[C] vector<lower=0>[D] sd_G; //sd of ind effects
  
  // season RE

  // // plot RE
  matrix[P, D] Z_nestbox;
  cholesky_factor_corr[D] L_nestbox;
  vector<lower=0>[D] sd_nestbox;

  // ZI coefficient
  real<lower=0, upper=1> theta;
  
  // cutpoints for ordinal model (equivalent to intercepts)
  ordered[14] cutpoint; 
}

model {
  //predicted values from reaction norms
  //growth
  vector[N] mu_growth =  Q_g * B_mq_g;
  
  //fecundity
  vector[M] mu_fecundity =  Q_f * B_mq_f;
  
  //recruitment
  vector[M] mu_recruitment =  Q_r * B_mq_r;
                       
  //correlations (expressed as canonical partial correlations)
  matrix[C, ncor] cpc_G = tanh(Q * B_cpcq);
  
  //variances
  matrix[C, D] sd_G = sqrt(exp(Q * B_vq));
  
  //scale univariate random effects
  
  //scale multivariate nestbox random effects
  matrix[P, D] mat_nestbox = Z_nestbox * diag_pre_multiply(sd_nestbox, L_nestbox)';


 //initialize mean linear predictors
  vector[N] mu_g = mu_growth[id_g_lm] + col(mat_nestbox, 1)[plot_id_g]; 
  vector[M] mu_f = mu_fecundity[id_f_lm] + col(mat_nestbox, 2)[plot_id_f]; 
  vector[M] mu_r = mu_recruitment[id_f_lm] + col(mat_nestbox, 3)[plot_id_f]; 

  //scale context-specific multivariate additive genetic effects
  matrix[cnt, D] mat_G;
  int pos = 1; //keep track of position 1:cnt
  for(c in 1:C){
      mat_G[pos:(pos+cn[c]-1)] = 
      LA[cmat[c,1:cn[c]],cmat[c,1:cn[c]]] * Z_G[pos:(pos+cn[c]-1)] * diag_pre_multiply(sd_G[c],lkj_to_chol_corr(cpc_G[c], D))';
      pos = pos + cn[c];   
  }
        
//add context-specific genetic effects to linear predictors
  for(n in 1:N){
  mu_g[n] += col(mat_G,1)[idc_g[n]];
  }
  
  for(m in 1:M){
  mu_f[m] += col(mat_G,2)[idc_f[m]]; 
  mu_r[m] += col(mat_G,3)[idc_f[m]]; 
  }        
       
       
//likelihood growth (gaussian)
  growth ~ normal(mu_g, sd_E);
//likelihood fecundity (poisson)
  //productivity ~ poisson_log(mu_f);
  productivity ~ ordered_logistic(mu_f, cutpoint);
//likelihood recruitment (zero inflated poisson)
   vector[M_zero] mu_r_zero = mu_r[id_zero];
   vector[M - M_zero] mu_r_nonzero = mu_r[id_nonzero];
  
   target += M_zero * log_sum_exp(log(theta),
                                  log1m(theta) + poisson_lpmf(0 | exp(mu_r_zero) ));
   target += M_nonzero * log1m(theta);
   target += poisson_lpmf(y_nonzero | exp(mu_r_nonzero));
  
  
//priors
  to_vector(B_mq_g) ~ normal(0,1);
  to_vector(B_mq_f) ~ normal(0,1);
  to_vector(B_mq_r) ~ normal(0,1);
  to_vector(B_cpcq) ~ normal(0,0.5);
  to_vector(B_vq) ~ normal(0,1);
  to_vector(Z_G) ~ std_normal();
  
  theta ~ beta(1,1);
  
  to_vector(Z_nestbox) ~ std_normal();
  sd_nestbox ~ exponential(2);
  L_nestbox ~ lkj_corr_cholesky(2);

  sd_E ~ exponential(2);
  
  // for(c in 1:C){
  // sd_G[c] ~ exponential(2);
  // }
}

generated quantities{
  vector[P_g] B_m_g; //mean RN parameters for X
  vector[P_f] B_m_f; //mean RN parameters for X
  vector[P_r] B_m_r; //mean RN parameters for X
  matrix[P_y,ncor] B_cpc; //partial correlation RN parameters for X
  matrix[P_y,D] B_v;

  B_m_g= R_inv_g * B_mq_g;
  B_m_f= R_inv_f * B_mq_f;
  B_m_r= R_inv_r * B_mq_r;

  for(d in 1:ncor){
    B_cpc[,d]= R_inv * B_cpcq[,d];
    }
    
  for(i in 1:D){
    B_v[,i]= R_inv * B_vq[,i];
    }  
    
  matrix[D, D] Sigma_plot = L_nestbox * L_nestbox';
  
    // Posterior predictive check for fecundity and offpsring mass
  vector[N] mu_growth_bis =  Q_g * B_mq_g;
  vector[M] mu_fecundity_bis =  Q_f * B_mq_f;
  vector[M] mu_recruitment_bis =  Q_r * B_mq_r;
  matrix[C, ncor] cpc_G_bis = tanh(Q * B_cpcq);
  matrix[C, D] sd_G_bis = sqrt(exp(Q * B_vq));
  
  matrix[P, D] mat_nestbox_bis = Z_nestbox * diag_pre_multiply(sd_nestbox, L_nestbox)';
  vector[N] mu_g_bis = mu_growth_bis[id_g_lm] + col(mat_nestbox_bis, 1)[plot_id_g];
  vector[M] mu_f_bis = mu_fecundity_bis[id_f_lm] + col(mat_nestbox_bis, 2)[plot_id_f];
  vector[M] mu_r_bis = mu_recruitment_bis[id_f_lm] + col(mat_nestbox_bis, 3)[plot_id_f];
  matrix[cnt, D] mat_G_bis;
  int pos = 1; //keep track of position 1:cnt
  for(c in 1:C){
      mat_G_bis[pos:(pos+cn[c]-1)] =
      LA[cmat[c,1:cn[c]],cmat[c,1:cn[c]]] * Z_G[pos:(pos+cn[c]-1)] * diag_pre_multiply(sd_G_bis[c],lkj_to_chol_corr(cpc_G_bis[c], D))';
      pos = pos + cn[c];
  }

  array[N] real y_rep_g;
  array[M] int y_rep_f;
  array[M] int y_rep_r;
  for(n in 1:N){
  mu_g_bis[n] += col(mat_G_bis,1)[idc_g[n]];
  }
  for(m in 1:M){
  mu_f_bis[m] += col(mat_G_bis,2)[idc_f[m]];
  mu_r_bis[m] += col(mat_G_bis,3)[idc_f[m]];
  }

  y_rep_g = normal_rng(mu_g_bis, sd_E);
  //y_rep_f = poisson_log_rng(mu_f_bis);
  for(m in 1:M){
  y_rep_f[m] = ordered_logistic_rng(mu_f_bis[m], cutpoint);
  y_rep_r[m] = bernoulli_logit_rng(theta) ? 0 : poisson_log_rng(mu_r_bis[m]);
  }

}

