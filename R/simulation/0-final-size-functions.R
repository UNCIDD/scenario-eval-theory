#### FUNCTIONS -----------------------------------------------------------------
#' function to implement full simulation of final size simulation experiment
#' 
#' @description
#' The full simulation includes (1) generation of final size 
#' scenario projections for n_models across n_locations under two vaccination 
#' coverage scenarios; (2) calculating the error various errors for each model
#' and the truth
#' 
#' @param n_locations integer, number of locations to simulate
#' @param n_models integer, number of models to simulate
#' @param vax_cov_S1 double, assumed vaccination coverage in scenario 1 (lower value)
#' @param vax_cov_S2 double, assumed vaccination coverage in scenario 2 (upper value)
#' @param true_vax_cov_lwr double, lower bound on true vaccination coverage values
#'                         if NA, assumed to be the same as vax_cov_S1
#' @param true_vax_cov_upr double, upper bound on true vaccination coverage values
#'                         if NA, assumed to be the same as vax_cov_S2
#' @param lengthout_scenarios number of scenario values to return results for (minimum = 2, for S1 and S2)
#' @param R0_lwr double, lower bound on uniform distribution from which location-specific R0 is drawn
#' @param R0_upr double, upper bound on uniform distribution from which location-specific R0 is drawn
#' @param cov_R0 SOMETHING ABOUT CORRELATIONS IN OBSERVATIONS
#' @param model_bias_R0_mean double, individual models are assigned a mean bias
#'                           on R0; this specifies the mean of the distribution
#'                           from which this model-specific bias is drawn
#' @param model_bias_R0_sd double, individual models are assigned a mean bias on
#'                         R0; this specifies the standard deviation of the distribution
#'                         from which this model-specific bias is drawn
#' @param model_bias_ind_sd double, for each model, the location-specific bias in 
#'                          R0 is drawn from a normal distribution with model bias
#'                          as the mean (drawn using params `model_bias_R0_mean` 
#'                          and `model_bias_R0_sd`); this specifies the standard
#'                          deviation of that distribution (i.e., how consistent
#'                          the bias of an individual model in estimating R0
#'                          across locations); defaults to 0, i.e., no variation 
#'                          in model bias across locations
#' @param final_size_method "analytical" to use final_size() function, "simulation"
#'                          to use deterministic SIR model with I^alpha
#' @param alpha_lwr double lower bound for uniform distribtuion to draw model-specific alpha
#' @param alpha_upr double upper bound for uniform distribution to draw model-specific alpha
#' @param alpha_sd double standard deviation to draw location-specific alpha (around model-specific alpha)
#' 
#' 
#' 
#' @details
#' The simulation proceeds in the following steps: 
#' 1. generate predictions of final size from each model under specified scenarios 
#'    and true values for each location
#' 2. calculate errors for each projection, across locations/models
#' 
#' To generate predictions of final size, we first draw a realized vaccination
#' coverage and a true R0 value for each location. We calculate the true final
#' epidemic size for the low and high vaccination scenarios, as well as the
#' realized vaccination scenario. If cov_R0 is NA, these values are drawn
#' independently from uniform distributions (i.e., for each location i,
#' true_vax_cov_i ~ U(`true_vax_cov_lwr`, `true_vax_cov_upr`) and 
#' R0_i ~ U(`R0_lwr`, `R0_upr`)). However, if a value is specified for `cov_R0`,
#' we  draw correlated values for realized vaccination coverage and true R0 from
#' a multi-variate normal distribution, with covariance `cov_R0`. 
#' 
#' We calculate the true final epidemic size for the low and high vaccination 
#' scenarios, as well as the realized vaccination scenario. There are two possible 
#' methods to calculate final size, `final_size_method = "analytical"` which
#' uses the `finalSize` package, or `final_size_method = "simulation"` which
#' solves a system of ODEs over 1.5 years.
#' 
#' Then, after true values for each location have been drawn, we draw a model-
#' specific bias in R0 estimates, where for model j, 
#' model_bias_j ~ N(`model_bias_R0_mean`, `model_bias_R0_sd`). This provides
#' control of whether models tend to over- or under-estimate R0 across
#' locations. The model estimated R0 in a given location is the true R0 for 
#' that location plus the model bias and the location specific bias. 
#' 
#' We also include the possibility of scaling the infection term 
#' I^alpha. Again, we draw a model-specific alpha, where for model j, 
#' model_alpha_j ~ U(`alpha_lwr`, `alpha_upr`) and the location-specific alpha is 
#' then alpha_i ~ N(model_alpha_j, alpha_sd). We set the true alpha to 
#' `mean(alpha_lwr, alpha_upr)`. 
#' 
#' The projected final size is again calculated based on model estimated R0 and alpha 
#' for both scenarios of interest and for the realized vaccination coverage value.
#' 
#' Once projections and true values are generated, we calculate the errors for 
#' each projection as projected final size minus true final size. The entire 
#' process is also performed for intermediate values between scenarios so the 
#' true error relationship for each location is available if desired. These are
#' recorded with scenario_id = "E".
#' 
#' @return list, including simulated values for each model/location/scenario 
#' (model_sims), true values for each location (true_sims), and errors for each 
#' model/location/scenario (errors)
full_sim <- function(
    n_locations, n_models, seed = 100, vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
    true_vax_cov_lwr = NA, true_vax_cov_upr = NA, lengthout_scenarios = 20,
    R0_lwr = 2, R0_upr = 3.5, 
    cov_R0 = NA, model_bias_R0_mean = 0, model_bias_R0_sd = 0.05, 
    model_bias_ind_sd = 0, final_size_method = "analytical", 
    alpha_upr = 1, alpha_lwr = 0.95, alpha_sd = 0.01){
  if(is.na(true_vax_cov_lwr)){true_vax_cov_lwr = vax_cov_S1}
  if(is.na(true_vax_cov_upr)){true_vax_cov_upr = vax_cov_S2}
  sims <- generate_final_size_preds(n_locations, n_models, seed, 
                                    true_vax_cov_lwr, true_vax_cov_upr, cov_R0,
                                    vax_cov_S1, vax_cov_S2, lengthout_scenarios, 
                                    R0_lwr, R0_upr, 
                                    model_bias_R0_mean, model_bias_R0_sd, 
                                    model_bias_ind_sd, final_size_method, 
                                    alpha_upr, alpha_lwr, alpha_sd)
  model_sims = sims %>% filter(model_id != "T") %>%
    dplyr::select(-location_R0)
  true_sims = sims %>% filter(model_id == "T") %>%
    rename(true_final_size = final_size) %>%
    dplyr::select(location_id, scenario_id, location_R0, vax_cov, true_final_size)
  errors <- calculate_errors(model_sims, true_sims)$errors
  return(list(model_sims = model_sims, true_sims = true_sims, errors = errors))
}

#### HELPERS -------------------------------------------------------------------

#' function to simulate final size predictions from each model and the truth
generate_final_size_preds <- function(n_locations, n_models, seed, 
                                      true_vax_cov_lwr, true_vax_cov_upr, cov_R0,
                                      vax_cov_S1, vax_cov_S2, lengthout_scenarios,
                                      R0_lwr, R0_upr, 
                                      model_bias_R0_mean, model_bias_R0_sd, 
                                      model_bias_ind_sd, final_size_method, 
                                      alpha_upr, alpha_lwr, alpha_sd){
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
  # bias in R0 estimate for each model and location
  # first get model bias
  model_bias_R0 = data.frame(
    model_id = c(paste0("M", 1:n_models), "T"),
    model_bias = c(rnorm(n_models, model_bias_R0_mean, model_bias_R0_sd), 0) # 0 for true model
  )
  # then combine into an R0 value for each model-location pair (if model_id == T, there is no bias)
  model_loc_R0 = expand.grid(model_id =  c(paste0("M", 1:n_models), "T"), 
                             location_id = 1:n_locations) %>% 
    left_join(model_bias_R0) %>% 
    left_join(data.frame(location_id = 1:n_locations, 
                         location_R0 = R0_T)) %>% 
    mutate(location_bias_sd = ifelse(model_id == "T", 0, model_bias_ind_sd), 
           R0 = location_R0 + rnorm(1, model_bias, location_bias_sd), .by = c("model_id", "location_id"))
  # add alpha variable (if needed)
  if(final_size_method == "simulation"){
    model_loc_R0 = model_loc_R0 %>%
      left_join(
        data.frame(model_id = c(paste0("M", 1:n_models), "T"),
                   model_alpha = c(runif(n_models, alpha_lwr, alpha_upr), mean(c(alpha_lwr, alpha_upr)))) # ensure the true alpha is in the middle
      ) %>%
      mutate(alpha = rnorm(1, model_alpha, alpha_sd), .by = c("model_id", "location_id"))
  }
  else{
    model_loc_R0$alpha = 1
  }
  # generate all possibilities
  sims <- expand.grid(model_id = c(paste0("M", 1:n_models), "T"),
                      location_id = 1:n_locations, 
                      scenario_id = c("S1", "T", "S2"))
  sims$vax_cov = with(sims, ifelse(scenario_id == "T", vax_cov_T[location_id], 
                                   ifelse(scenario_id == "S1", vax_cov_S1, vax_cov_S2)))
  # generate extra sims along scenario axis to see true relationship
  if(lengthout_scenarios > 2){
    sims_full_relationship <- expand.grid(model_id =  c(paste0("M", 1:n_models), "T"), 
                                          location_id = 1:n_locations, 
                                          vax_cov = seq(vax_cov_S1, vax_cov_S2, length.out = lengthout_scenarios)) %>%
      mutate(scenario_id = "E")
    sims <- bind_rows(sims, sims_full_relationship)
  }
  sims <- sims %>%
    # add model bias
    left_join(model_loc_R0[, c("model_id", "location_id", "location_R0", "R0", "alpha")])
  sims$final_size = NA
  # final size variables
  susc_immunised <- cbind(1,0)
  colnames(susc_immunised) <- c("novax", "vax")
  for(i in 1:nrow(sims)){
    p_susc_immunised <- cbind(
      susceptible = 1 - sims$vax_cov[i],
      immunised = sims$vax_cov[i]
    )
    if(final_size_method == "analytical"){
      fs <- final_size(
        r0 = sims$R0[i],
        contact_matrix = matrix(1.0)/1E3,
        demography_vector = 1E3,
        susceptibility = susc_immunised,
        p_susceptibility = p_susc_immunised, 
        control = list(iterations = 1e8)
      )
      # browser()
      sims$final_size[i] <- fs[1,4]
    }
    else if(final_size_method == "simulation"){
      params = c(mu = 0, N = 1e3, R0 = unname(unlist(sims$R0[i])), gamma = 365/10, 
                 alpha = unname(unlist(sims$alpha[i])))
      params["beta"] = unname(unlist(params["R0"] * (params["gamma"] + params["mu"])))
      inits = c(S = 1 - unname(unlist(sims$vax_cov[i])),
                I = 1e-3,
                R = unname(unlist(sims$vax_cov[i])) - 1e-3
      )*params["N"]
      fs <- as.data.frame(ode(y = inits, times = seq(0, 1.5, 1/52), func = ode_sir, parms = params))
      sims$final_size[i] <- (inits["S"] - unname(unlist(fs %>% filter(time == max(time)) %>% pull(S))))/(inits["S"])
    }
  }
  return(sims)
}

#' function to calculate errors for each model and the truth
calculate_errors <- function(model_sims, true_sims){
  errors = model_sims %>% 
    left_join(true_sims, by = join_by(location_id, scenario_id, vax_cov)) %>%
    mutate(error = final_size - true_final_size, 
           relerror = (final_size - true_final_size)/true_final_size)
  return(list(model_sims = model_sims, true_sims = true_sims, errors = errors))
}

#' SIR differential equations
ode_sir = function(t, y, parameters) {
  with(as.list(c(y, parameters)), {
    # Define equations
    dS = mu * N - beta * S * I^alpha/N - mu * S
    dI = beta * S * I^alpha/N - gamma * I - mu * I
    dR = gamma * I - mu * R 
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  })
}

