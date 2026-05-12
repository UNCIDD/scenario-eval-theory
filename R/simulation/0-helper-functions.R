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
                         xval = as.vector(xvals)), by = join_by(xval_id)) %>% 
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


#' Calculate coverage of a given distribtuion
#' 
#' @param 
calculate_coverage <- function(quant_fits, error_df, vax_cov_S1, vax_cov_S2, summarize_by = "all"){
  cov_quant <- quant_fits %>%
    filter(vax_cov %in% c(vax_cov_S1, vax_cov_S2)) %>%
    mutate(scenario_id = ifelse(vax_cov == vax_cov_S1, "S1", "S2")) %>%
    filter(quantile != 0.5) %>%
    mutate(alpha = round(ifelse(quantile < 0.5, 1-2*quantile, 1-2*(1-quantile)),3), 
           bound = ifelse(quantile < 0.5, "lwr", "upr")) %>%
    reshape2::dcast(model_id + vax_cov + scenario_id + alpha ~ bound, value.var = "value") %>%
    dplyr::select(model_id, scenario_id, alpha, lwr, upr) %>%
    left_join(
      error_df %>%
        filter(scenario_id %in% c("S1", "S2")) %>%
        mutate(obs = error) %>%
        dplyr::select(model_id, location_id, scenario_id, obs),
      relationship = "many-to-many", by = join_by(model_id, scenario_id)
    ) %>%
    mutate(cov = ifelse(obs <= upr & obs >= lwr, 1, 0))
  if(summarize_by == "all"){
    cov_quant <- cov_quant  %>%
      summarize(cov = sum(cov)/n(), .by = c("scenario_id", "alpha")) 
  }
  else if(summarize_by == "model"){
    cov_quant <- cov_quant  %>%
      summarize(cov = sum(cov)/n(), .by = c("scenario_id", "model_id", "alpha")) 
  }
  return(cov_quant)
}

