data {
  int<lower=1> n_sites; //number of sites in study
  int<lower=1> n_years; //number of years of study  
  int<lower=1> n_obs; //total number of observations  
  array[n_obs] int<lower=0> dets; // Observations ordered by site -> replicate -> year
  array[n_obs] int<lower=0> site_yr_idx; // array of site-years indices to keep track
  int<lower=1> n_occ_covs;                // number of covariates for occupancy
  matrix[n_sites, n_occ_covs] xmat; //(scaled) covariate matrix for occupancy
  int<lower=1> n_det_covs;                // number of covariates for detection
  matrix[n_obs, n_det_covs] vmat;         // (scaled) covariate matrix for detection

  //Priors
  real phi_mu; real phi_sigma;
  real gamma_mu; real gamma_sigma;
  real alpha0_mu; real alpha0_sigma; 
  array[n_occ_covs] real alphas_mu; array[n_occ_covs] real alphas_sigma;
  real beta0_mu; real beta0_sigma; 
  array[n_det_covs] real betas_mu; array[n_det_covs] real betas_sigma;
}

parameters {
  // Set up model parameters
  vector[n_years - 1] logit_phi; // survival
  vector[n_years - 1] logit_gamma; // recruitment

  // Initial occupancy parameters
  real alpha0;                  // Baseline intercept for occupancy (Year 1)
  vector[n_occ_covs] alphas;    // Coefficients for environmental covariates (xmat)
  
  // Nightly detection parameters
  real beta0;                   // Intercept for true-positive detection
  vector[n_det_covs] betas;     // Coefficients for nightly survey covariates (vmat)
}

transformed parameters {
  // Occupancy Model (Year 1)
  vector[n_sites] logit_psi1 = alpha0 + xmat * alphas;
  vector[n_sites] psi1 = inv_logit(logit_psi1);

  vector<lower = 0, upper = 1> [n_years - 1] phi = inv_logit(logit_phi);
  vector<lower = 0, upper = 1> [n_years - 1] gamma = inv_logit(logit_gamma);
}

model {
    // Priors
    logit_phi ~ normal(phi_mu, phi_sigma);
    logit_gamma ~ normal(gamma_mu, gamma_sigma);
    alpha0 ~ normal(alpha0_mu, alpha0_sigma);
    alphas ~ normal(alphas_mu, alphas_sigma);
    beta0 ~ normal(beta0_mu, beta0_sigma);
    betas ~ normal(betas_mu, betas_sigma);

    // Detection Model
    vector[n_obs] logit_p = beta0 + vmat * betas;

    // Set up TPM
    array[n_years - 1] matrix[2, 2] log_tpm_t;
    matrix[2,2] log_tpm_temp;
    for (t in 1:n_years - 1) {
        log_tpm_temp[1, 1] = log(phi[t]);
        log_tpm_temp[1, 2] = log1m(phi[t]);
        log_tpm_temp[2, 1] = log(gamma[t]);
        log_tpm_temp[2, 2] = log1m(gamma[t]);

        log_tpm_t[t] = log_tpm_temp';
    }
    
    // Initialize Pointers
    int pos_obs = 1;
    int sy_idx = 1; 

    for (i in 1:n_sites) {

        // Detection Model Log-likelihood Year 1
        vector[2] lp;
        real ll_occ_y1 = 0.0;
        real ll_unocc_y1 = 0.0;
        
        while (pos_obs <= n_obs && site_yr_idx[pos_obs] == sy_idx) {
            ll_occ_y1 += bernoulli_logit_lpmf(dets[pos_obs] | logit_p[pos_obs]);
            if (dets[pos_obs] == 1) {
                ll_unocc_y1 += negative_infinity(); //No False Positives Constraint
            }
            pos_obs += 1;
        }

        // Compute log-prob for year 1
        lp[1] = log_inv_logit(logit_psi1[i]) + ll_occ_y1;
        lp[2] = log1m_inv_logit(logit_psi1[i]) + ll_unocc_y1;
        
        sy_idx += 1;

        // Repeat for remaining years
        for (t in 2:n_years) {
        real ll_occ = 0.0;
        real ll_unocc = 0.0;
        
        while (pos_obs <= n_obs && site_yr_idx[pos_obs] == sy_idx) {
            ll_occ += bernoulli_logit_lpmf(dets[pos_obs] | logit_p[pos_obs]);
            if (dets[pos_obs] == 1) {
            ll_unocc += negative_infinity();
            }
            pos_obs += 1;
        }

        vector[2] lp_p1;
        lp_p1[1] = log_sum_exp(to_vector(log_tpm_t[t - 1, 1]) + lp) + ll_occ;
        lp_p1[2] = log_sum_exp(to_vector(log_tpm_t[t - 1, 2]) + lp) + ll_unocc;
        lp = lp_p1;
        
        sy_idx += 1;
        }

        target += log_sum_exp(lp);
    }
}

generated quantities {
    // Variables declared here will be saved in the final output
    matrix<lower = 0, upper = 1>[n_sites, n_years] psi;
    array[n_obs] int <lower = 0, upper = 1> y_rep;

    // Compute expected psi in year 1
    psi[,1] = psi1;

    // Compute expected psi in year t
    for (i in 1:n_sites){
        for (t in 2:n_years){ 
            psi[i,t] = psi[i,t-1] * phi[t-1] + (1 - psi[i,t-1]) * gamma[t-1];
        }
    }

    // --- LOCAL BLOCK: Variables declared here are discarded ---
    {
        // Declare z_sim locally so it isn't saved
        array[n_sites, n_years] int z_sim;
        
        // Detection probs declared locally
        vector[n_obs] logit_p = beta0 + vmat * betas;
        vector[n_obs] p = inv_logit(logit_p);
        
        // Simulate the latent states (z_sim)
        for (i in 1:n_sites) {
            // Year 1 initial state
            z_sim[i, 1] = bernoulli_rng(psi1[i]);

            // Years 2 to T transitions
            for (t in 2:n_years) {
              real prob_occ = z_sim[i, t - 1] * phi[t - 1] + (1 - z_sim[i, t - 1]) * gamma[t - 1];
              z_sim[i, t] = bernoulli_rng(prob_occ);
            }
        }
        
        // Generate posterior predictive observations (y_rep) using sampled z
        for (n in 1:n_obs) {
            int sy = site_yr_idx[n];
            
            // Map linear site-year index back to (site, year)
            int s = (sy - 1) %/% n_years + 1;
            int y_idx = (sy - 1) % n_years + 1;

            // Detection probability depends on occupancy: if z_sim == 0, prob is 0
            real eff_p = z_sim[s, y_idx] * p[n];
            y_rep[n] = bernoulli_rng(eff_p);
        }
    }
}

