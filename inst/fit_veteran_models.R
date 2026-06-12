# fit_veteran_models.R
# Analysis script to investigate differences between approaches
# Saves results to R/sysdata.rda (preserves other objects like getting_started_results)

library(brms)
library(survival)
library(shrinkr)
library(distributional)
library(tidybayes)
library(dplyr)
library(posterior)

# Set seed for reproducibility
set.seed(1104)

# Load data
data(veteran, package = "survival")
cell_types <- c("squamous", "smallcell", "adeno", "large")

cat(strrep("=", 72), "\n")
cat("VETERAN DATASET ANALYSIS - INVESTIGATING DIFFERENCES\n")
cat(strrep("=", 72), "\n\n")

# ============================================================================
# APPROACH 1: Two-Stage (brms + shrinkr)
# ============================================================================

cat("APPROACH 1: Two-stage (brms + shrinkr)\n")
cat(strrep("-", 72), "\n\n")

cat("Step 1a: Fitting brms Cox model (uninformative priors)...\n")
cat("Model: time | cens(1 - status) ~ trt:celltype + karno + age\n\n")

brms_uninformative <- brm(
  time | cens(1 - status) ~ trt:celltype + karno + age,
  data = veteran,
  family = cox(),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = 0.95),
  refresh = 0,
  silent = 2
)

# Extract summary
brms_uninformative_summary <- capture.output(print(summary(brms_uninformative)))

cat("Step 1b: Extracting posteriors...\n")
brms_posteriors <- brms_uninformative %>%
  spread_draws(`b_trt:celltypesquamous`, 
               `b_trt:celltypesmallcell`,
               `b_trt:celltypeadeno`,
               `b_trt:celltypelarge`) %>%
  select(-c(.chain, .iteration, .draw)) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = "celltype",
    values_to = "value"
  ) %>%
  mutate(celltype = gsub("b_trt:celltype", "", celltype)) %>%
  group_by(celltype) %>%
  summarise(draws = list(matrix(value, ncol = 1)), .groups = "drop") %>%
  tibble::deframe()

# Reorder to match cell_types vector
brms_posteriors <- brms_posteriors[cell_types]

cat("Checking posterior order:\n")
cat("Expected order:", paste(cell_types, collapse = ", "), "\n")
cat("Actual order:", paste(names(brms_posteriors), collapse = ", "), "\n\n")

# Check posterior correlations (KEY DIAGNOSTIC)
cat("Step 1c: Checking posterior correlations...\n")
stage1_mat <- do.call(cbind, brms_posteriors)
colnames(stage1_mat) <- cell_types
stage1_cov <- cov(stage1_mat)
stage1_cor <- cor(stage1_mat)

cat("Posterior correlations between cell types:\n")
print(round(stage1_cor, 3))
cat("\nMean off-diagonal correlation:", 
    round(mean(stage1_cor[upper.tri(stage1_cor)]), 3), "\n\n")

cat("Step 1d: Fitting mixture...\n")
mix_brms <- fit_mixture(
  samples = brms_posteriors,
  K_max = 3,
  verbose = TRUE
)

# ============================================================================
# APPROACH 2: Full brms hierarchical (single stage)
# ============================================================================

cat("\n\nAPPROACH 2: Full brms hierarchical (single stage)\n")
cat(strrep("-", 72), "\n\n")

cat("Fitting brms Cox model with hierarchical structure...\n")
cat("Model: time | cens(1 - status) ~ trt + (trt | celltype) + karno + age\n")
cat("Using Half-Normal(0, 0.5) prior on sd to match shrinkr\n\n")

brms_hierarchical <- brm(
  time | cens(1 - status) ~ trt + (0 + trt | celltype) + karno + age,
  data = veteran,
  family = cox(),
  prior = c(
    prior(normal(0, 1), class = b, coef = "trt"),
    prior(normal(0, 0.5), class = sd, group = celltype, lb = 0)  # Half-Normal
  ),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = 0.95),
  refresh = 0,
  silent = 2
)

# Extract summaries
brms_hierarchical_summary <- capture.output(print(summary(brms_hierarchical)))

# Extract cell-specific effects (b_trt + r_celltype)
brms_hier_effects <- brms_hierarchical %>%
  spread_draws(b_trt, r_celltype[celltype, term]) %>%
  filter(term == "trt") %>%
  mutate(effect = b_trt + r_celltype) %>%
  group_by(celltype) %>%
  summarise(
    mean = mean(effect),
    sd = sd(effect),
    q2.5 = quantile(effect, 0.025),
    q50 = quantile(effect, 0.50),
    q97.5 = quantile(effect, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    hr_mean = exp(mean),
    hr_lower = exp(q2.5),
    hr_upper = exp(q97.5)
  )

# Extract tau (sd_celltype__trt)
brms_tau <- brms_hierarchical %>%
  spread_draws(sd_celltype__trt) %>%
  summarise(
    mean = mean(sd_celltype__trt),
    sd = sd(sd_celltype__trt),
    q2.5 = quantile(sd_celltype__trt, 0.025),
    q97.5 = quantile(sd_celltype__trt, 0.975)
  )

# Extract mu (b_trt)
brms_mu <- brms_hierarchical %>%
  spread_draws(b_trt) %>%
  summarise(
    mean = mean(b_trt),
    sd = sd(b_trt),
    q2.5 = quantile(b_trt, 0.025),
    q97.5 = quantile(b_trt, 0.975)
  )

cat("brms hierarchical tau (sd_celltype__trt):\n")
print(brms_tau)
cat("\nbrms hierarchical mu (b_trt):\n")
print(brms_mu)

# ============================================================================
# APPROACH 3: Frequentist Cox + shrinkr
# ============================================================================

cat("\n\nAPPROACH 3: Frequentist Cox + shrinkr\n")
cat(strrep("-", 72), "\n\n")

cat("Fitting frequentist Cox model...\n")
cat("Model: Surv(time, status) ~ trt:celltype + karno + age\n\n")

cox_model <- coxph(
  Surv(time, status) ~ trt:celltype + karno + age,
  data = veteran
)

# Extract summary
cox_summary <- summary(cox_model)

# Extract treatment effects
coefs <- coef(cox_model)
vcov_mat <- vcov(cox_model)
trt_terms <- grep("trt:celltype", names(coefs))
trt_effects <- coefs[trt_terms]
trt_vcov <- vcov_mat[trt_terms, trt_terms]
names(trt_effects) <- gsub("trt:celltype", "", names(trt_effects))
rownames(trt_vcov) <- colnames(trt_vcov) <- names(trt_effects)

# Reorder to match cell_types
trt_effects <- trt_effects[cell_types]
trt_vcov <- trt_vcov[cell_types, cell_types]

cat("Frequentist estimates (log HR):\n")
print(trt_effects)
cat("\nFrequentist standard errors:\n")
print(sqrt(diag(trt_vcov)))

# ============================================================================
# SENSITIVITY ANALYSIS: 5 different tau priors
# ============================================================================

cat("\n\nSENSITIVITY ANALYSIS\n")
cat(strrep("-", 72), "\n\n")

# Define 5 half-normal priors
prior_specs <- list(
  very_strong = list(
    name = "Very Strong",
    label = "Half-Normal(0, 0.1)",
    scale = 0.1,
    priors = list(
      mu = dist_normal(0, 1),
      tau = dist_truncated(dist_normal(0, 0.1), lower = 0)
    )
  ),
  strong = list(
    name = "Strong",
    label = "Half-Normal(0, 0.25)",
    scale = 0.25,
    priors = list(
      mu = dist_normal(0, 1),
      tau = dist_truncated(dist_normal(0, 0.25), lower = 0)
    )
  ),
  moderate = list(
    name = "Moderate",
    label = "Half-Normal(0, 0.5)",
    scale = 0.5,
    priors = list(
      mu = dist_normal(0, 1),
      tau = dist_truncated(dist_normal(0, 0.5), lower = 0)
    )
  ),
  weak = list(
    name = "Weak",
    label = "Half-Normal(0, 1.0)",
    scale = 1.0,
    priors = list(
      mu = dist_normal(0, 1),
      tau = dist_truncated(dist_normal(0, 1.0), lower = 0)
    )
  ),
  very_weak = list(
    name = "Very Weak",
    label = "Half-Normal(0, 2.0)",
    scale = 2.0,
    priors = list(
      mu = dist_normal(0, 1),
      tau = dist_truncated(dist_normal(0, 2.0), lower = 0)
    )
  )
)

# Fit all sensitivity models
sensitivity_summaries <- list()

for (spec_name in names(prior_specs)) {
  spec <- prior_specs[[spec_name]]
  cat(sprintf("\nFitting %s borrowing (%s)...\n", spec$name, spec$label))
  
  # Fit with brms posteriors
  fit_brms <- shrink(
    mixture = mix_brms,
    hierarchical_priors = spec$priors,
    chains = 4,
    iter = 4000,
    warmup = 1000,
    cores = 4,
    verbose = FALSE
  )
  
  sensitivity_summaries[[paste0(spec_name, "_brms")]] <- list(
    theta_summary = summarise_theta(fit_brms, group_names = cell_types),
    mu_tau_summary = summarise_mu_tau(fit_brms),
    print_output = capture.output(print(fit_brms))
  )
  
  # Fit with frequentist estimates
  fit_freq <- shrink(
    mle = trt_effects,
    var_matrix = trt_vcov,
    hierarchical_priors = spec$priors,
    chains = 4,
    iter = 4000,
    warmup = 1000,
    cores = 4,
    verbose = FALSE
  )
  
  sensitivity_summaries[[paste0(spec_name, "_freq")]] <- list(
    theta_summary = summarise_theta(fit_freq, group_names = cell_types),
    mu_tau_summary = summarise_mu_tau(fit_freq),
    print_output = capture.output(print(fit_freq))
  )
}

# ============================================================================
# KEY COMPARISONS
# ============================================================================

cat("\n\n", strrep("=", 72), "\n")
cat("KEY COMPARISONS AND DIAGNOSTICS\n")
cat(strrep("=", 72), "\n\n")

# Compare theta estimates
comparison <- data.frame(
  celltype = cell_types,
  brms_hierarchical = brms_hier_effects$mean[match(cell_types, brms_hier_effects$celltype)],
  brms_shrinkr = sensitivity_summaries$moderate_brms$theta_summary$mean,
  freq_shrinkr = sensitivity_summaries$moderate_freq$theta_summary$mean
)

comparison$diff_brms_shrinkr <- comparison$brms_hierarchical - comparison$brms_shrinkr
comparison$diff_brms_freq <- comparison$brms_hierarchical - comparison$freq_shrinkr

cat("THETA ESTIMATES (log HR scale):\n")
print(comparison)

cat("\nMaximum absolute difference:\n")
cat("brms hierarchical vs brms+shrinkr:", 
    round(max(abs(comparison$diff_brms_shrinkr)), 4), "\n")
cat("brms hierarchical vs freq+shrinkr:", 
    round(max(abs(comparison$diff_brms_freq)), 4), "\n")

# Compare tau estimates
cat("\n\nTAU ESTIMATES:\n")
cat("\nbrms hierarchical (sd_celltype__trt):\n")
print(brms_tau)

cat("\nshrinkr (from brms stage 1, moderate prior):\n")
print(sensitivity_summaries$moderate_brms$mu_tau_summary %>% filter(parameter == "tau"))

cat("\nshrinkr (from freq stage 1, moderate prior):\n")
print(sensitivity_summaries$moderate_freq$mu_tau_summary %>% filter(parameter == "tau"))

# Compare mu estimates
cat("\n\nMU ESTIMATES:\n")
cat("\nbrms hierarchical (b_trt):\n")
print(brms_mu)

cat("\nshrinkr (from brms stage 1, moderate prior):\n")
print(sensitivity_summaries$moderate_brms$mu_tau_summary %>% filter(parameter == "mu"))

cat("\nshrinkr (from freq stage 1, moderate prior):\n")
print(sensitivity_summaries$moderate_freq$mu_tau_summary %>% filter(parameter == "mu"))

# Stage 1 estimates comparison
cat("\n\nSTAGE 1 ESTIMATES (before shrinkage):\n")
brms_stage1_means <- sapply(brms_posteriors, mean)
cat("\nbrms stage 1 posterior means:\n")
print(round(brms_stage1_means, 4))

cat("\nFrequentist MLEs:\n")
print(round(trt_effects, 4))

cat("\nDifference (brms - freq):\n")
print(round(brms_stage1_means - trt_effects, 4))

# Correlation check
cat("\n\nPOSTERIOR CORRELATIONS (explains why approaches may differ):\n")
print(round(stage1_cor, 3))

# ============================================================================
# SAVE TO R/sysdata.rda — save both objects together so neither overwrites the other
# ============================================================================

cat("\n\nSaving results to R/sysdata.rda...\n")

# Strip any lingering brms / tidybayes class attributes that R CMD check
# would otherwise flag as a hidden namespace dependency on brms. Casting
# tibble-like objects through tibble::as_tibble(as.data.frame(.)) gives
# us a plain tibble without inherited classes.
.scrub <- function(x) {
  if (is.data.frame(x)) tibble::as_tibble(as.data.frame(x))
  else if (is.list(x))  lapply(x, .scrub)
  else                  x
}
brms_posteriors    <- lapply(brms_posteriors, function(m) `attributes<-`(m, list(dim = dim(m))))
brms_hier_effects  <- .scrub(brms_hier_effects)
brms_tau           <- .scrub(brms_tau)
brms_mu            <- .scrub(brms_mu)

veteran_analysis <- list(
  cell_types = cell_types,
  brms_uninformative_summary = brms_uninformative_summary,
  brms_posteriors = brms_posteriors,
  stage1_cor = stage1_cor,
  stage1_cov = stage1_cov,
  mix_brms = mix_brms,
  brms_hierarchical_summary = brms_hierarchical_summary,
  brms_hier_effects = brms_hier_effects,
  brms_tau = brms_tau,
  brms_mu = brms_mu,
  cox_summary = cox_summary,
  trt_effects = trt_effects,
  trt_vcov = trt_vcov,
  prior_specs = prior_specs,
  sensitivity_summaries = sensitivity_summaries,
  comparison = comparison
)

sysdata_path <- "R/sysdata.rda"

# Load getting_started_results if it exists, otherwise save veteran_analysis alone
if (file.exists(sysdata_path)) {
  e <- new.env()
  load(sysdata_path, envir = e)
  if (exists("getting_started_results", envir = e)) {
    getting_started_results <- get("getting_started_results", envir = e)
    cat("   - Loaded existing getting_started_results from sysdata.rda\n")
    usethis::use_data(
      veteran_analysis,
      getting_started_results,
      internal  = TRUE,
      overwrite = TRUE,
      compress  = "xz"
    )
  } else {
    cat("   - No getting_started_results found, saving veteran_analysis only\n")
    usethis::use_data(
      veteran_analysis,
      internal  = TRUE,
      overwrite = TRUE,
      compress  = "xz"
    )
  }
} else {
  cat("   - No existing sysdata.rda, saving veteran_analysis only\n")
  usethis::use_data(
    veteran_analysis,
    internal  = TRUE,
    overwrite = TRUE,
    compress  = "xz"
  )
}

final_size <- file.size(sysdata_path)
cat("   - Saved. Final size:", format(final_size, units = "MB"), "\n")

# Verify both objects are present
e_check <- new.env()
load(sysdata_path, envir = e_check)
cat("   - Objects in sysdata.rda:", paste(ls(e_check), collapse = ", "), "\n")

cat("\n", strrep("=", 72), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 72), "\n")