#' Simulate prediction intervals from GAM fit
#' 
#' adapted from: https://www.mail-archive.com/r-help@r-project.org/msg132608.html
#' example with normal: https://mikl.dk/post/2019-prediction-intervals-for-gam/
#' here assuming normal distribution
#' 
#' @param mod GAM model fit
#' @param xvals vector of x values for which to return prediction intervals
#' @param invfn2 function, additional inverse function to be applied if predictions
#' were transformed before fitting GAM, otherwise identity use function
get_gam_PIs <- function(mod, xvals, invfn2 = function(x){return(x)}){
  beta <- coef(mod) # beta
  Vb <- vcov(mod) # V
  # simulate beta vectors 
  reps <- 10000 # num_beta_vecs
  nb <- length(beta)
  br <- t(chol(Vb)) %*% matrix(rnorm(reps*nb), nb, reps) + beta # beta_sims <- beta + t(Cv) %*% matrix(nus, length_beta, num_beta_vecs) 
  # replicates to linear predictors
  Xp <- predict(mod, newdata = xvals, type = "lpmatrix") # covar_sim (using grid of xvals instead of random samples)
  lp <- Xp%*%br # linpred_sim = covar_sim %*% beta_sims
  invfun <- family(mod)$linkinv # invlink
  fv <- invfun(lp) # exp_val_sim
  yr <- matrix(rnorm(fv*0, mean = fv, sd = sqrt(summary(mod)$scale)), # y_sim
               nrow(fv), ncol(fv)) 
  # transform and summarize
  ret <- reshape2::melt(yr, keep.rownames = TRUE) %>%
    rename(xval_id = Var1, sim = Var2) %>% 
    mutate(value_inv = invfn2(value)) %>% 
    reframe(quantile = quantiles, 
            value = quantile(value_inv, quantiles), .by = c("xval_id")) %>% 
    left_join(data.frame(xval_id = 1:length(unlist(xvals)), 
                         xval = xvals)) %>% 
    dplyr::select(-xval_id)
  return(ret)
}

#' Get samples from a distribution
#' 
#' @param quantile vector of quantiles
#' @param value vector of values
#' @param n_samp integer with number of samples
#' @param seed
get_samps <- function(quantile, value, n_samps = 1e4, seed = 1002){
  set.seed(seed)
  approx(quantile, value, runif(n_samps), yleft = min(value), yright = max(value))$y
}
