// stage2_shrinkage.stan
// ============================================================================
// Two-Stage Bayesian Hierarchical Shrinkage Model
// ============================================================================
//
// This Stan program implements the Stage 2 hierarchical model for shrinkr.
// It takes as input a Gaussian mixture approximation of Stage 1 posteriors
// and applies hierarchical shrinkage with flexible prior specifications.
//
// Model:
//   theta_g ~ Normal(mu, tau)         [hierarchical model]
//   theta_g ~ q(theta_g | Stage 1)    [mixture approximation likelihood]
//
// Supported priors:
//   - Built-in: Normal, Student-t, Lognormal, Inverse-Gamma, Gamma,
//               Exponential, Uniform, Cauchy
//   - Mixture priors (e.g., spike-and-slab)
//   - Truncated priors on mu
//
// Parameters:
//   mu:  Global mean effect across groups
//   tau: Between-group heterogeneity (standard deviation)
//   theta[g]: Group-specific effects (g = 1, ..., G)
//
// Stan array syntax: uses the `array[N] type name` form introduced in
// Stan 2.26 (stanc3). Requires rstan >= 2.26 / StanHeaders >= 2.26.
//
// ============================================================================

functions {
  // Evaluate built-in priors (precompiled)
  real builtin_tau_prior_lpdf(real tau, int prior_code, array[] real params) {
    if (prior_code == 1) {
      // Normal: params[1]=mu, params[2]=sigma
      return normal_lpdf(tau | params[1], params[2]);

    } else if (prior_code == 2) {
      // Student-t: params[1]=mu, params[2]=sigma, params[3]=df
      return student_t_lpdf(tau | params[3], params[1], params[2]);

    } else if (prior_code == 3) {
      // Lognormal: params[4]=mu, params[5]=sigma
      return lognormal_lpdf(tau | params[4], params[5]);

    } else if (prior_code == 4) {
      // Inverse-Gamma: params[1]=shape, params[2]=scale
      return inv_gamma_lpdf(tau | params[1], params[2]);

    } else if (prior_code == 5) {
      // Gamma: params[1]=shape, params[2]=rate
      return gamma_lpdf(tau | params[1], params[2]);

    } else if (prior_code == 6) {
      // Exponential: params[1]=rate
      return exponential_lpdf(tau | params[1]);

    } else if (prior_code == 7) {
      // Uniform: handled by bounds (tau_lb, tau_ub)
      return 0;  // uniform density is constant within bounds

    } else if (prior_code == 8) {
      // Cauchy: params[1]=location, params[2]=scale
      return cauchy_lpdf(tau | params[1], params[2]);
    }

    return 0; // Should not reach
  }

  // Evaluate mixture prior (works for both mu and tau)
  real mixture_normal_prior_lpdf(real x, int n_components, vector weights,
                                 vector locs, vector scales) {
    vector[n_components] lps;
    for (j in 1:n_components) {
      lps[j] = log(weights[j]) + normal_lpdf(x | locs[j], scales[j]);
    }
    return log_sum_exp(lps);
  }
}

data {
  // Dimension
  int<lower=1> G;                     // number of groups
  int<lower=1> K;                     // mixture components for q(theta)

  // Mixture q(theta) as DATA
  simplex[K] w;                       // weights
  array[K] vector[G] m;               // means
  array[K] cholesky_factor_cov[G] L;  // Cholesky factors

  // ===== MU PRIOR SPECIFICATION =====
  int<lower=1, upper=3> mu_prior_type;  // 1=Normal, 2=Student-t, 3=Mixture
  real mu_loc;
  real<lower=0> mu_scale;
  real<lower=0> mu_df;                  // For Student-t

  // For mixture prior on mu
  int<lower=1> mu_n_components;
  simplex[mu_n_components] mu_mix_weights;
  vector[mu_n_components] mu_mix_locs;
  vector<lower=0>[mu_n_components] mu_mix_scales;

  // Truncation bounds for mu
  real mu_lb;
  real mu_ub;
  int<lower=0, upper=1> mu_is_truncated;

  // ===== TAU PRIOR SPECIFICATION =====
  // Prior type: 1-8 = built-in, 9 = mixture, 10 = custom
  int<lower=1, upper=10> tau_prior_type;

  // Parameters for built-in priors
  // [1]=loc/shape, [2]=scale/rate, [3]=df, [4]=lognormal_mu, [5]=lognormal_sigma, [6]=lower_bound
  array[6] real tau_params;

  // For mixture prior on tau
  int<lower=1> tau_n_components;
  simplex[tau_n_components] tau_mix_weights;
  vector[tau_n_components] tau_mix_locs;
  vector<lower=0>[tau_n_components] tau_mix_scales;

  // Bounds for tau
  real<lower=0> tau_lb;
  real<lower=tau_lb> tau_ub;

  // Custom prior specification
  int<lower=0> custom_prior_type;  // 0=none, 1=simple expression, 2=complex
  array[10] real custom_params;    // Parameters for custom prior

  // Parameterization
  int<lower=0, upper=1> centered;
}

parameters {
  real<lower=(mu_is_truncated ? mu_lb : negative_infinity()),
       upper=(mu_is_truncated ? mu_ub : positive_infinity())> mu;
  real<lower=tau_lb, upper=tau_ub> tau;

  // Group effects: only one set is "active" depending on parameterization
  vector[G] theta_c;             // used if centered==1
  vector[G] z;                   // used if centered==0
}

transformed parameters {
  vector[G] theta;
  if (centered == 1) {
    theta = theta_c;
  } else {
    theta = mu + tau * z;
  }
}

model {
  // ===== PRIOR ON MU =====
  if (mu_prior_type == 1) {
    mu ~ normal(mu_loc, mu_scale);
  } else if (mu_prior_type == 2) {
    mu ~ student_t(mu_df, mu_loc, mu_scale);
  } else if (mu_prior_type == 3) {
    target += mixture_normal_prior_lpdf(mu | mu_n_components, mu_mix_weights,
                                        mu_mix_locs, mu_mix_scales);
  }

  // ===== PRIOR ON TAU =====
  if (tau_prior_type <= 8) {
    target += builtin_tau_prior_lpdf(tau | tau_prior_type, tau_params);
  } else if (tau_prior_type == 9) {
    target += mixture_normal_prior_lpdf(tau | tau_n_components, tau_mix_weights,
                                        tau_mix_locs, tau_mix_scales);
  } else if (tau_prior_type == 10) {
    if (custom_prior_type == 1) {
      target += -custom_params[1] * tau - log(2 * custom_params[1]);
    } else if (custom_prior_type == 2) {
      target += -custom_params[1] * log(tau) - custom_params[2] / tau;
    }
  }

  // ===== HIERARCHICAL MODEL =====
  // Both parameterizations get a prior so the inactive one is constrained.
  // Without this, the inactive vector has an improper uniform on (-inf, inf)
  // which causes Stan to waste time exploring unbounded space.
  if (centered == 1) {
    theta_c ~ normal(mu, tau);
    z ~ std_normal();             // constrains unused z
  } else {
    z ~ std_normal();
    theta_c ~ std_normal();       // constrains unused theta_c
  }

  // ===== STAGE-2 LIKELIHOOD =====
  if (K == 1) {
    target += multi_normal_cholesky_lpdf(theta | m[1], L[1]);
  } else {
    vector[K] lps;
    for (k in 1:K) {
      lps[k] = log(w[k]) + multi_normal_cholesky_lpdf(theta | m[k], L[k]);
    }
    target += log_sum_exp(lps);
  }
}

generated quantities {
  real<lower=0> tau_squared = square(tau);
}
