## M: set size

##---------------
##  Slots Model  
##---------------

## k: capacity (number of slots)
## g: guess change 

## m: prob. item in memory = min(1, k/M)
## p(Hit) = m + (1 - m)*g
## p(FA) = (1 - m)*g

# set up custom brms family:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

slots <- custom_family(
  "slots", dpars = c("mu", "g"),
  links = c("log", "logit"),
  lb = c(0, 0), ub = c(NA, 1),
  type = "int", vars = c("vint1[n]", "vint2[n]")
)

stan_funs_slots <- "
  real slots_lpmf(int y, real mu, real g, int type, int setsize) {
    real m;
    real p;
    m = min({1, mu/setsize});
    if (type == 0) {
      p = (1 - m)*g;
    } else if (type == 1) {
      p = m + (1 - m)*g;
    }
    return bernoulli_lpmf(y | p);
  }
"
stanvars_slots <- stanvar(scode = stan_funs_slots, block = "functions")

calc_posterior_predictions_slots <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  g <- brms::get_dpar(prep, "g", i = i)
  OUTLEN <- length(mu)
  type <- prep$data$vint1[i]
  setsize <- prep$data$vint2[i]
  p <- vector("numeric", OUTLEN)
  
  m <- pmin(1, mu/setsize)

  if (type == 0) {
    p = (1 - m)*g  
  } else if (type == 1) {
    p = m + (1 - m)*g
  }
  return(p)
}

log_lik_slots <- function(i, prep) {
  p <- calc_posterior_predictions_slots(i = i, prep = prep)
  out <- p
  out[prep$data$Y[i] == 0] <- 1 - p
  log(out)
}

posterior_epred_slots <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs)
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_slots(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

posterior_predict_slots <- function(i, prep, ...) {
  p <- calc_posterior_predictions_slots(i = i, prep = prep)
  rbinom(length(p), 1, p)
}
