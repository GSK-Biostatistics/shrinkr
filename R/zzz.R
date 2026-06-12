.onLoad <- function(libname, pkgname) {
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    try(Rcpp::loadModule("stan_fit4stage2_shrinkage_mod", what = TRUE), silent = TRUE)
  }
}

.onAttach <- function(...) {
  stanv <- rstan::stan_version()
  packageStartupMessage(paste("shrinkr (Version ", utils::packageVersion("shrinkr"),
                               ") using rstan (Version ", stanv, ")", sep = ""))
}
