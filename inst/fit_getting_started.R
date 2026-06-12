# Generate pre-computed results for getting started vignette
# Only saves Stage 1 outputs (slow Stan fit, mixture, prior predictive)
# shrink() runs live in the vignette so nothing from it is stored here
# Location: inst/generate_getting_started_data.R

library(shrinkr)
library(rstan)
library(distributional)
library(dplyr)

set.seed(1104)

cat(strrep("=", 78), "\n")
cat("Generating pre-computed data for getting_started vignette\n")
cat(strrep("=", 78), "\n\n")

# ============================================================================
# 1. Generate trial data
# ============================================================================
cat("1. Generating trial data...\n")

true_effects <- c(0.45, 0.72, 0.38, 0.55, 0.61)
regions      <- c("North", "South", "East", "West", "Central")
n_per_region <- c(100, 80, 120, 90, 70)

trial_data <- lapply(seq_along(regions), function(i) {
  n <- n_per_region[i]
  data.frame(
    region    = regions[i],
    treatment = rep(c(0, 1), each = n/2),
    outcome   = c(
      rnorm(n/2, mean = 0, sd = 1),
      rnorm(n/2, mean = true_effects[i], sd = 1)
    )
  )
}) %>% bind_rows()

cat("   - Total observations:", nrow(trial_data), "\n\n")

# ============================================================================
# 2. Fit Stan model (the slow part)
# ============================================================================
cat("2. Fitting Stan model...\n")

stan_code <- "
data {
  int<lower=0> N;
  int<lower=1> G;
  vector[N] y;
  vector[N] treatment;
  array[N] int<lower=1,upper=G> region;
}
parameters {
  vector[G] beta_region;
  real<lower=0> sigma;
}
model {
  sigma ~ normal(0, 2);
  for (n in 1:N) {
    y[n] ~ normal(treatment[n] * beta_region[region[n]], sigma);
  }
}
"

region_indices <- as.integer(factor(trial_data$region, levels = regions))

stan_fit <- stan(
  model_code = stan_code,
  data = list(
    N = nrow(trial_data), G = length(regions),
    y = trial_data$outcome, treatment = trial_data$treatment,
    region = region_indices
  ),
  chains = 4, iter = 2000, warmup = 1000, refresh = 0, seed = 123
)

beta_draws   <- extract(stan_fit)$beta_region
samples_list <- lapply(seq_along(regions), function(i) beta_draws[, i])
names(samples_list) <- regions
samples <- lapply(samples_list, function(x) matrix(x, ncol = 1))

# Discard Stan fit object immediately - it carries colorspace dependency
rm(stan_fit)
gc(verbose = FALSE)

cat("   - Extracted", nrow(samples[[1]]), "draws per region\n\n")

# ============================================================================
# 3. Fit mixture
# ============================================================================
cat("3. Fitting mixture approximation...\n")

mix <- fit_mixture(samples = samples, K_max = 3, verbose = FALSE)

cat("   - Selected K =", mix$K, "components\n\n")

# ============================================================================
# 4. Generate prior predictive samples
# ============================================================================
cat("4. Generating prior predictive samples...\n")

hierarchical_priors <- list(
  mu  = dist_normal(0, 1),
  tau = dist_truncated(dist_student_t(3, 0, 0.5), lower = 0)
)

prior_pred <- sample_prior_predictive(
  hierarchical_priors = hierarchical_priors,
  n_groups = 5,
  n_draws  = 1000
)

cat("   - Generated", nrow(prior_pred$theta), "prior draws\n\n")

# ============================================================================
# 5. Build getting_started_results
# ============================================================================
cat("5. Building getting_started_results...\n")

getting_started_results <- list(
  trial_data = trial_data,
  samples    = samples,
  mix        = mix,
  prior_pred = prior_pred
)

cat("\n   Object sizes:\n")
cat("   - trial_data:", format(object.size(trial_data),  units = "KB"), "\n")
cat("   - samples:   ", format(object.size(samples),     units = "MB"), "\n")
cat("   - mix:       ", format(object.size(mix),         units = "KB"), "\n")
cat("   - prior_pred:", format(object.size(prior_pred),  units = "KB"), "\n")
cat("   - TOTAL:     ", format(object.size(getting_started_results), units = "MB"), "\n\n")

# ============================================================================
# 6. Save both objects together so neither overwrites the other
# ============================================================================
cat("6. Saving to R/sysdata.rda...\n")

sysdata_path <- "R/sysdata.rda"

# Load veteran_analysis if it exists, otherwise use NULL as placeholder
if (file.exists(sysdata_path)) {
  e <- new.env()
  load(sysdata_path, envir = e)
  if (exists("veteran_analysis", envir = e)) {
    veteran_analysis <- get("veteran_analysis", envir = e)
    cat("   - Loaded existing veteran_analysis from sysdata.rda\n")
    usethis::use_data(
      getting_started_results,
      veteran_analysis,
      internal  = TRUE,
      overwrite = TRUE,
      compress  = "xz"
    )
  } else {
    cat("   - No veteran_analysis found, saving getting_started_results only\n")
    usethis::use_data(
      getting_started_results,
      internal  = TRUE,
      overwrite = TRUE,
      compress  = "xz"
    )
  }
} else {
  cat("   - No existing sysdata.rda, saving getting_started_results only\n")
  usethis::use_data(
    getting_started_results,
    internal  = TRUE,
    overwrite = TRUE,
    compress  = "xz"
  )
}

cat("   - Saved. Final size:", format(file.size(sysdata_path), units = "MB"), "\n")

# Verify both objects are present
e_check <- new.env()
load(sysdata_path, envir = e_check)
cat("   - Objects in sysdata.rda:", paste(ls(e_check), collapse = ", "), "\n")

cat("\n", strrep("=", 78), "\n")
cat("Done! Re-run this script whenever Stage 1 data needs to be refreshed.\n")
cat(strrep("=", 78), "\n")