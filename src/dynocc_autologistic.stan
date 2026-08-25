data {
  // Dimensions
  int<lower=1> n_sites; 
  int<lower=1> n_years;   
  int<lower=1> n_obs;   
  
  // Observation data
  array[n_obs] int<lower=0> dets; 
  array[n_sites, n_years] int<lower=0> start_idx; 
  array[n_sites, n_years] int<lower=0> end_idx;   
  
  // Covariates
  int<lower=1> n_occ_covs;                
  matrix[n_sites, n_occ_covs] xmat; 
  int<lower=1> n_det_covs;                
  matrix[n_obs, n_det_covs] vmat;         
  
  // Prior hyperparameters 
  real a_mu; 
  real<lower=0> a_sigma;
  
  real b_mu; 
  real<lower=0> b_sigma;
  
  real alpha0_mu; 
  real<lower=0> alpha0_sigma; 
  
  vector[n_occ_covs] alphas_mu; 
  vector<lower=0>[n_occ_covs] alphas_sigma;
  
  real beta0_mu; 
  real<lower=0> beta0_sigma; 
  
  vector[n_det_covs] betas_mu; 
  vector<lower=0>[n_det_covs] betas_sigma;
}

parameters {
  // Transition parameters (time-varying)
  vector[n_years - 1] a; 
  vector[n_years - 1] b; 
  
  // Occupancy parameters (spatial)
  real alpha0;                  
  vector[n_occ_covs] alphas;    
  
  // Detection parameters
  real beta0;                   
  vector[n_det_covs] betas;     
}

transformed parameters {
  // Pre-compute spatial effects 
  vector[n_sites] psi_spatial = alpha0 + xmat * alphas;
  vector[n_sites] trans_spatial = xmat * alphas;
}

model {
  // Priors supplied dynamically via data block
  a ~ normal(a_mu, a_sigma);
  b ~ normal(b_mu, b_sigma);
  alpha0 ~ normal(alpha0_mu, alpha0_sigma);
  alphas ~ normal(alphas_mu, alphas_sigma);
  beta0 ~ normal(beta0_mu, beta0_sigma);
  betas ~ normal(betas_mu, betas_sigma);

  // Nightly detection logit
  vector[n_obs] logit_p = beta0 + vmat * betas;

  // Forward Algorithm
  for (i in 1:n_sites) {
    vector[2] lp;
    real ll_occ_y1 = 0.0;
    real ll_unocc_y1 = 0.0;
    
    int s1 = start_idx[i, 1];
    int e1 = end_idx[i, 1];
    
    // Log-likelihood for Year 1
    if (s1 > 0) {
      for (n in s1:e1) {
        ll_occ_y1 += bernoulli_logit_lpmf(dets[n] | logit_p[n]);
        if (dets[n] == 1) {
          ll_unocc_y1 = negative_infinity(); 
        }
      }
    }

    // Initialize the forward variable
    lp[1] = log_inv_logit(psi_spatial[i]) + ll_occ_y1;
    lp[2] = log1m_inv_logit(psi_spatial[i]) + ll_unocc_y1;

    // Log-likelihood for remaining years
    for (t in 2:n_years) {
      real ll_occ_t = 0.0;
      real ll_unocc_t = 0.0;
      
      int st = start_idx[i, t];
      int et = end_idx[i, t];
      
      // Log-likelihood in year T
      if (st > 0) {
        for (n in st:et) {
          ll_occ_t += bernoulli_logit_lpmf(dets[n] | logit_p[n]);
          if (dets[n] == 1) {
            ll_unocc_t = negative_infinity();
          }
        }
      }

      // Site-and-time-specific transition logits
      real logit_gamma = a[t-1] + trans_spatial[i];
      real logit_phi = a[t-1] + b[t-1] + trans_spatial[i];

      vector[2] lp_next;
      
      // Marginalize over the latent states
      lp_next[1] = log_sum_exp(
          lp[1] + log_inv_logit(logit_phi),     
          lp[2] + log_inv_logit(logit_gamma)    
      ) + ll_occ_t;
      
      lp_next[2] = log_sum_exp(
          lp[1] + log1m_inv_logit(logit_phi),   
          lp[2] + log1m_inv_logit(logit_gamma)  
      ) + ll_unocc_t;

      lp = lp_next;
    }
    
    // Increment log-likelihood for the site
    target += log_sum_exp(lp);
  }
}

generated quantities {
  matrix<lower=0, upper=1>[n_sites, n_years] psi;
  vector<lower=0, upper=1>[n_years - 1] phi;
  vector<lower=0, upper=1>[n_years - 1] gamma;
  array[n_obs] int<lower=0, upper=1> y_rep;

  // Exported variables for residual diagnostics
  array[n_sites, n_years] int z_sim;
  matrix[n_sites, n_years] occ_res;
  vector[n_obs] p;
  vector[n_obs] det_res;

  // Derived regional average transition parameters
  // Evaluated at average environmental conditions 
  for (t in 1:(n_years - 1)) {
    phi[t] = inv_logit(a[t] + b[t]);
    gamma[t] = inv_logit(a[t]);
  }

  //Compute site-specific expected occupancy (psi)
  for (i in 1:n_sites) {
    psi[i, 1] = inv_logit(psi_spatial[i]);
    
    for (t in 2:n_years) {
      real phi_it = inv_logit(a[t-1] + b[t-1] + trans_spatial[i]);
      real gamma_it = inv_logit(a[t-1] + trans_spatial[i]);
      
      psi[i, t] = psi[i, t-1] * phi_it + (1 - psi[i, t-1]) * gamma_it;
    }
  }

  //Reconstruct detection probability for all observations
  p = inv_logit(beta0 + vmat * betas); 

  //Simulate latent states and calculate occupancy residuals
  for (i in 1:n_sites) {
    // Year 1
    z_sim[i, 1] = bernoulli_rng(psi[i, 1]);
    occ_res[i, 1] = z_sim[i, 1] - psi[i, 1]; 
    
    // Years 2 to T
    for (t in 2:n_years) {
      real prob_occ = z_sim[i, t-1] * inv_logit(a[t-1] + b[t-1] + trans_spatial[i]) + 
                      (1 - z_sim[i, t-1]) * inv_logit(a[t-1] + trans_spatial[i]);
                      
      z_sim[i, t] = bernoulli_rng(prob_occ);
      occ_res[i, t] = z_sim[i, t] - psi[i, t]; 
    }
  }

  //Generate y_rep
  for (i in 1:n_sites) {
    for (t in 1:n_years) {
      int st = start_idx[i, t];
      int et = end_idx[i, t];
      
      if (st > 0) {
        for (n in st:et) {
          y_rep[n] = bernoulli_rng(z_sim[i, t] * p[n]);
          det_res[n] = y_rep[n] - p[n];
        }
      }
    }
  }
}