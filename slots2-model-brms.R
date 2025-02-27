## M: set size

##---------------------------
##  Slots + Attention Model  
##---------------------------

## k: capacity (number of slots)
## g: guess change 

## m: prob. item in memory = min(1, k/M)
## p(Hit) = a*m + a*(1 - m)*g + (1-a)*g
## p(FA) = a*(1 - m)*g + (1-a)*g

# set up custom brms family:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

slots2 <- custom_family(
  "slots2", dpars = c("mu", "g", "a"),
  links = c("log", "logit", "logit"),
  lb = c(0, 0, 0), ub = c(NA, 1, 1),
  type = "int", vars = c("vint1[n]", "vint2[n]")
)

stan_funs_slots2 <- "
  real slots2_lpmf(int y, real mu, real g, real a, int type, int setsize) {
    real m;
    real p;
    m = min({1, mu/setsize});
    if (type == 0) {
      p = a*(1 - m)*g + (1-a)*g;
    } else if (type == 1) {
      p = a*m + a*(1 - m)*g + (1-a)*g;
    }
    return bernoulli_lpmf(y | p);
  }
"
stanvars_slots2 <- stanvar(scode = stan_funs_slots2, block = "functions")

calc_posterior_predictions_slots2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  g <- brms::get_dpar(prep, "g", i = i)
  a <- brms::get_dpar(prep, "a", i = i)
  OUTLEN <- length(mu)
  type <- prep$data$vint1[i]
  setsize <- prep$data$vint2[i]
  p <- vector("numeric", OUTLEN)
  
  m <- pmin(1, mu/setsize)

  if (type == 0) {
    p = a*(1 - m)*g + (1-a)*g
  } else if (type == 1) {
    p = a*m + a*(1 - m)*g + (1-a)*g
  }
  return(p)
}

log_lik_slots2 <- function(i, prep) {
  p <- calc_posterior_predictions_slots2(i = i, prep = prep)
  out <- p
  out[prep$data$Y[i] == 0] <- 1 - p
  log(out)
}

posterior_epred_slots2 <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs)
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_slots2(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

posterior_predict_slots2 <- function(i, prep, ...) {
  p <- calc_posterior_predictions_slots2(i = i, prep = prep)
  rbinom(length(p), 1, p)
}
