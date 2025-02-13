#### FUNCTIONS -----------------------------------------------------------------
#' function to implement full simulation of final size simulation experiment
#' 
#' @description
#' The full simulation includes (1) generation of final size 
#' scenario projections for n_models across n_locations under two vaccination 
#' coverage scenarios; (2) estimating the error distribution for each model, and
#' (3) calculating the coverage of the estimated error distribution against the
#' true error distribution for each model
#' 
#' @param n_locations integer, number of locations to simulate
#' @param n_models integer, number of models to simulate
#' @param vax_cov_S1 double, assumed vaccination coverage in scenario 1 (lower value)
#' @param vax_cov_S2 double, assumed vaccination coverage in scenario 2 (upper value)
#' @param true_vax_cov_lwr double, lower bound on true vaccination coverage values
#'                         if NA, assumed to be the same as vax_cov_S1
#' @param true_vax_cov_upr double, upper bound on true vaccination coverage values
#'                         if NA, assumed to be the same as vax_cov_S2
#' @param R0_lwr double, lower bound on uniform distribution from which location-specific R0 is drawn
#' @param R0_upr double, upper bound on uniform distribution from which location-specific R0 is drawn
#' @param cov_R0 SOMETHING ABOUT CORRELATIONS IN OBSERVATIONS
#' @param model_bias_R0_mean ADD HERE
#' @param model_bias_R0_sd ADD HERE
#' @param fit_outcomes logical, TRUE to estimate error distribution and calculate coverage
#' 
#' @details
#' The simulation proceeds in the following steps: 
#' 1. generate predictions of final size from each model under specified scenarios 
#'    and true values for each location
#' 2. calculate errors for each projection, across locations/models
#' 3. (optional) estimate error distribution 
#' 
#' To generate predictions of final size, we first draw a realized vaccination
#' coverage and a true R0 value for each location. Using the `finalSize` package, 
#' we calculate the true final epidemic size for the low and high vaccination 
#' scenarios, as well as the realized vaccination scenario. If cov_R0 is NA, 
#' these values are drawn independently from uniform distributions (i.e., for 
#' each location i, true_vax_cov_i ~ U(`true_vax_cov_lwr`, `true_vax_cov_upr`) and 
#' R0_i ~ U(`R0_lwr`, `R0_upr`)). However, if a value is specified for `cov_R0`, we 
#' draw correlated values for realized vaccination coverage and true R0 from
#' a multi-variate normal distribution, with covariance `cov_R0`. 
#' 
#' Then, after true values for each location have been drawn, we draw a model-
#' sepcific bias in R0 estimates, where for model j, 
#' model_bias_j ~ N(`model_bias_R0_mean`, `model_bias_R0_sd`). This provides
#' control of whether models tend to over- or under-estimate R0 for a given 
#' location. Then, the model estimated R0 in a given location is the true R0 for 
#' that location plus the model bias. The projected final size is again 
#' calculated using the `finalSize` package based on model estimated R0 for 
#' both scenarios of interest and for the realized vaccination coverage value.
#' 
#' Once projections and true values are generated, we calculate the errors for 
#' each projection as projected final size minus true final size. The entire 
#' process is also performed for intermediate values between scenarios so the 
#' true error relationship for each location is available if desired. These are
#' recorded with scenario_id = "E".
#' 
#' TO ADD/UPDATE: DETAILS ABOUT ESTIMATING ERROR DISTRIBUTION
#' 
#' @return list, including simulated values for each model/location/scenario 
#' (model_sims), true values for each location(true_sims), and errors for each 
#' model/location/scenario (errors); if error distribution is estimated, a 
#' data.frame containing estimated quantiles for each model and scnenario (error_df)
#' and the corresponding coverage (coverage) are also included in the list
full_sim <- function(
    n_locations, n_models, seed = 100, vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
    true_vax_cov_lwr = NA, true_vax_cov_upr = NA, R0_lwr = 2, R0_upr = 3.5, 
    cov_R0 = NA, model_bias_R0_mean = 0, model_bias_R0_sd = 0.05, 
    fit_outcomes = TRUE){
  if(is.na(true_vax_cov_lwr)){true_vax_cov_lwr = vax_cov_S1}
  if(is.na(true_vax_cov_upr)){true_vax_cov_upr = vax_cov_S2}
  sims <- generate_final_size_preds(n_locations, n_models, seed, 
                                    true_vax_cov_lwr, true_vax_cov_upr, cov_R0,
                                    vax_cov_S1, vax_cov_S2, R0_lwr, R0_upr, 
                                    model_bias_R0_mean, model_bias_R0_sd)
  model_sims = sims %>% filter(model_id != "T")
  true_sims = sims %>% filter(model_id == "T") %>%
    rename(true_final_size = final_size) %>%
    dplyr::select(location_id, scenario_id, vax_cov, true_final_size)
  errors <- calculate_errors(model_sims, true_sims)$errors
  if(fit_outcomes){
    quant_reg_est <- estimate_w_gamlss(vax_cov_S1 = vax_cov_S1, 
                                       vax_cov_S2 = vax_cov_S2, errors_df = errors, 
                                       n_models = n_models)
    cov <- calculate_coverage(quant_reg_est, errors, vax_cov_S1, vax_cov_S2)
    return(list(model_sims = model_sims, true_sims = true_sims, errors = errors, 
                quant_reg = quant_reg_est, coverage = cov))
  }
  else{
    return(list(model_sims = model_sims, true_sims = true_sims, errors = errors))
  }
}

#' 
generate_final_size_preds <- function(n_locations, n_models, seed, 
                                      true_vax_cov_lwr, true_vax_cov_upr, cov_R0,
                                      vax_cov_S1, vax_cov_S2, R0_lwr, R0_upr, 
                                      model_bias_R0_mean, model_bias_R0_sd){
  set.seed(seed)
  if(is.na(cov_R0)){
    # true vaccination coverage for each location
    vax_cov_T = runif(n_locations, true_vax_cov_lwr, true_vax_cov_upr)
    # true R0 for each location
    R0_T = runif(n_locations, R0_lwr, R0_upr)
  }
  else{
    sigma<-rbind(c(1,cov_R0), c(cov_R0,1))
    # create the mean vector
    mu<-c(0,0)
    # generate the multivariate normal distribution
    df<-as.data.frame(mvrnorm(n=n_locations, mu=mu, Sigma=sigma))
    true_vax_cov_lwr2 = runif(1, true_vax_cov_lwr, true_vax_cov_lwr*1.05)
    true_vax_cov_upr2 = runif(1, true_vax_cov_upr*0.95, true_vax_cov_upr)
    R0_lwr2 = runif(1, R0_lwr, R0_lwr*1.05)
    R0_upr2 = runif(1, R0_upr*0.95, R0_upr)
    vax_cov_T = (df$V1 - min(df$V1))/(max(df$V1) - min(df$V1)) * (true_vax_cov_upr2 - true_vax_cov_lwr2) + true_vax_cov_lwr2
    R0_T = (df$V2 - min(df$V2))/(max(df$V2) - min(df$V2)) * (R0_upr2 - R0_lwr2) + R0_lwr2
  }
  # bias in R0 estimate for each location
  model_bias_R0 = c(rnorm(n_models, model_bias_R0_mean, model_bias_R0_sd), 0) # 0 for true model
  names(model_bias_R0) = c(paste0("M", 1:n_models), "T")
  # generate all possibilities
  sims <- expand.grid(model_id = c(paste0("M", 1:n_models), "T"),
                      location_id = 1:n_locations, 
                      scenario_id = c("S1", "T", "S2"))
  sims$vax_cov = with(sims, ifelse(scenario_id == "T", vax_cov_T[location_id], 
                                   ifelse(scenario_id == "S1", vax_cov_S1, vax_cov_S2)))
  # generate extra sims along scenario axis to see true relationship
  sims_full_relationship <- expand.grid(model_id =  c(paste0("M", 1:n_models), "T"), 
                                        location_id = 1:n_locations, 
                                        vax_cov = seq(vax_cov_S1, vax_cov_S2, length.out = 20)) %>%
    mutate(scenario_id = "E")
  sims <- bind_rows(sims, sims_full_relationship)
  # add model bias
  sims$R0 = R0_T[sims$location_id] + model_bias_R0[sims$model_id]
  sims$final_size = NA
  # final size variables
  susc_immunised <- cbind(1,0)
  colnames(susc_immunised) <- c("novax", "vax")
  for(i in 1:nrow(sims)){
    p_susc_immunised <- cbind(
      susceptible = 1 - sims$vax_cov[i],
      immunised = sims$vax_cov[i]
    )
    fs <- final_size(
      r0 = sims$R0[i],
      contact_matrix = matrix(1.0)/1E3,
      demography_vector = 1E3,
      susceptibility = susc_immunised,
      p_susceptibility = p_susc_immunised, 
      control = list(iterations = 1e8)
    )
    sims$final_size[i] <- fs[1,4]
  }
  return(sims)
}

calculate_errors <- function(model_sims, true_sims){
  errors = model_sims %>% 
    left_join(true_sims, by = join_by(location_id, scenario_id, vax_cov)) %>%
    mutate(error = final_size - true_final_size, 
           relerror = (final_size - true_final_size)/true_final_size)
  return(list(model_sims = model_sims, true_sims = true_sims, errors = errors))
}

estimate_w_gamlss <- function(vax_cov_S1, vax_cov_S2, errors_df, n_models,
                                        quantiles = c(0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99)){
  new_vax_cov = data.frame(vax_cov = seq(vax_cov_S1, vax_cov_S2, length.out = 100))
  quant_fits <- vector("list", n_models)
  for(i in 1:n_models){
    dat_filt = errors_df %>%
      filter(model_id == paste0("M", i), scenario_id == "T") 
    gam_fit <- try({
      gamlss(error ~ cs(vax_cov),
             sigma.formula = ~ cs(vax_cov), data = dat_filt,  family = NO, trace=FALSE)#, family = GA)
    }, silent = TRUE)
    ### quantile regression alternative
    # quant_fits[[i]] = as.data.frame()
    # qr_fit = rq(error ~ vax_cov, tau = quantiles, data = dat_filt)
    # quant_fits[[i]] = as.data.frame(predict(qr_fit, newdata = new_vax_cov)) %>%
    #   mutate(vax_id = 1:n()) %>%
    #   reshape2::melt(c("vax_id")) %>%
    #   mutate(vax_cov = new_vax_cov[vax_id, "vax_cov"], 
    #          variable = as.double(gsub("tau= ", "", variable))) %>%
    #   rename(quantile = variable)
    if(any(class(gam_fit) == "try-error")){
      gam_fit <- try({
        gamlss(error ~ pb(vax_cov), method=CG(), #method = mixed(2,10),
               sigma.formula = ~ pb(vax_cov), data = dat_filt,  family = NO)#, family = GA) family = SN1
      }, silent = TRUE)
      if(any(class(gam_fit) == "try-error")){
        print("error")
        quant_fits[[i]] = NA
        next
      }
    }
    g_pls = centiles.pred(gam_fit, type = "centiles", 
                          cent = quantiles*100,
                          xvalues = new_vax_cov$vax_cov, xname = "vax_cov", 
                          data = dat_filt)
    quant_fits[[i]] = reshape2::melt(g_pls, c("x")) %>%
      mutate(quantile = quantiles[variable]) %>%
      rename(vax_cov = x) %>%
      select(-variable)
  }
  quant_fits <- bind_rows(quant_fits, .id = "model_id") %>%
    mutate(model_id = paste0("M", model_id)) 
  return(quant_fits)
}

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
get_gam_PIs <- function(mod, xvals, invfn2){
  # get estimates and covariance matrix
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

