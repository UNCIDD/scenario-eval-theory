### simulation experiment to test scenario evaluation theory
### assume projections of final outbreak size across varying vaccination rates

#### SETUP ---------------------------------------------------------------------
library(finalsize)
library(dplyr)
library(reshape2)
library(ggplot2)
library(quantreg)
library(gamlss)
library(bestNormalize)
library(MASS)
source("./R/final-size-functions.R")

scenario_labs = c("low vax scenario", "high vax scenario")
names(scenario_labs) = c("S1", "S2")

#### GENERATE SIMULATIONS TO TEST ----------------------------------------------
reps = 1
set.seed(10)
seed_id = sample(1:1E6, reps)
n_models = 10

# quantiles of interest for summary error distribution
quantiles = c(0.001, 0.01, 0.025, seq(0.05, 0.95, 0.05), 0.975, 0.99, 0.999)
# and related alpha values 
alphas = sort(unique(sapply(quantiles, function(i){round(ifelse(i < 0.5, 1-2*i, 1-2*(1-i)),3)})))
# alphas = alphas[-which(alphas == 0)]

# new vaccination coverage
new_vax_cov = data.frame(vax_cov = seq(0.3, 0.5, length.out = 100))

# run simulation - using 500 locations for now
t_large <- full_sim(n_locations = 500, n_models = n_models,
                  vax_cov_S1 = 0.3, vax_cov_S2 = 0.5,
                  R0_lwr = 2, R0_upr = 3,
                  model_bias_ind_sd = 0.05,
                  seed = seed_id, fit_outcomes = FALSE)

# note: scenario_id = T returns only observed errors (not true errors to test against)
error_df = t_large$errors %>% filter(scenario_id == "T") 

#### TEST NORMALITY AND SOME TRANSFORMATIONS -----------------------------------
dat_transform <- error_df %>% 
  mutate(boxcox = sign(error) * bestNormalize::boxcox(abs(error))$x.t,
         yj = yeojohnson(error)$x.t, 
         ordnorm = orderNorm(error)$x.t, 
         cuberoot = sign(error)*abs(error)^(1/3), 
         log = sign(error)*log(abs(error)), 
         adjlog = log(1 + error - min(error)),
         sinh = asinh(error),
         .by = "model_id") %>% 
  dplyr::select(-relerror) %>%
  melt(c("model_id", "location_id", "scenario_id", "vax_cov", "R0", "final_size", "true_final_size"))

ggplot(data = dat_transform, aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(vars(paste(model_id, variable)), scales = "free", nrow = n_models) + 
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

ggplot(data = dat_transform, aes(sample = value, color = variable)) + 
  stat_qq(shape = 21) + 
  stat_qq_line() + 
  facet_wrap(vars(paste(model_id, variable)), scales = "free", nrow = n_models) + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

# box cox looks good, but we have to deal with the negative numbers

#### FUNCTIONS TO FIT MIXTURE DISTRIBUTION -------------------------------------
# proposal: fit positive and negative errors separately, using box-cox transform
# on x and -x, respectively
# then combine into a single distribution
# issues: 
#    1. cases where there are a very small number of negative or positive values
#       (this becomes more common for small n)
#    2. how to combine into a single distribution
#       (currently just taking a weighted average based on the proportion of
#       errors in the data that are negative or positive, but could imagine
#       estimating this parameter instead)
#    3. sometimes prediction intervals are out of bounds for box-cox inverse
#       transform (currently setting it to lower_bound + 1e-5...)

boxcox_transform <- function(y, lambda) {
  if (lambda == 0) log(y) else (y^lambda - 1) / lambda
} 

inverse_boxcox <- function(y_trans, lambda) {
  if (lambda == 0) exp(y_trans) else (lambda * y_trans + 1)^(1/lambda)
}

fit_one_boxcox <- function(error_df, new_vax_cov, alphas, neg_flag = FALSE){
  if(neg_flag){
    multiplier = -1
  }
  else{multiplier = 1}
  # find optimal transform
  bx = boxcox(multiplier*error ~ vax_cov, data = error_df, plotit = FALSE)
  lambda = with(bx, x[which.max(y)])
  # fit linear model
  lm = lm(boxcox_transform(multiplier*error, lambda) ~ vax_cov, data = error_df)
  # get prediction intervals
  pred_nominal <- list()
  for(j in 1:length(alphas)){
    pred <- predict(lm, newdata = new_vax_cov, level = alphas[j], interval = "prediction")
    if(neg_flag){# when negative, -upr becomes lower bound, and vice versa
      lwr_col = which(colnames(pred) == "upr");
      upr_col = which(colnames(pred) == "lwr")
      } 
    else{
      lwr_col = which(colnames(pred) == "lwr")
      upr_col = which(colnames(pred) == "upr")
    } 
    # check for possible errors in box-cox inverse transform
    lwr_pred = pred[, lwr_col]
    upr_pred = pred[, upr_col]
    if(lambda > 0){
      if(any(lwr_pred < -1/lambda) | any(upr_pred < -1/lambda)){
        print(paste("correcting: alpha =", alphas[j]))
      }
      lwr_pred[which(lwr_pred < -1/lambda)] = -1/lambda + 1e-5
      upr_pred[which(upr_pred < -1/lambda)] = -1/lambda + 1e-5
    }
    else{
      if(any(lwr_pred > -1/lambda) | any(upr_pred > -1/lambda)){
        print(paste("correcting: alpha =", alphas[j]))
      }
      lwr_pred[which(lwr_pred > -1/lambda)] = -1/lambda + 1e-5
      upr_pred[which(upr_pred > -1/lambda)] = -1/lambda + 1e-5
    }
    pred_nominal[[j]] <- data.frame(
        vax_cov = new_vax_cov,
        fit = multiplier*inverse_boxcox(pred[, "fit"], lambda),
        lwr = multiplier*inverse_boxcox(lwr_pred, lambda),
        upr = multiplier*inverse_boxcox(upr_pred, lambda), 
        alpha = alphas[j]
      )
  }
  pred_nominal = bind_rows(pred_nominal) %>% 
    reshape2::melt(c("vax_cov", "fit", "alpha")) %>%
    mutate(quantile = ifelse(variable == "lwr", (1-alpha)/2, 1-(1-alpha)/2))
  return(pred_nominal)
}

# fit_thresh: # obs neg or posis less than fit_thresh, ignore completely?
# to do: add plot flag if interested
fit_mixture_w_boxcox <- function(model_errors, new_vax_cov, alphas,
                                 plot_flag = FALSE, plot_title = NA,
                                 n_draws, fit_thresh){ 
    neg_errors <- model_errors %>% filter(error < 0)
    pos_errors <- model_errors %>% filter(error > 0)
    prop_neg = nrow(neg_errors)/nrow(model_errors)
    if(nrow(neg_errors) < fit_thresh){ # only fit positives
      print("fitting only positives")
      pred_nominal = fit_one_boxcox(pos_errors, new_vax_cov, alphas)
      pred_nominal = pred_nominal 
      r <- pred_nominal %>% 
        filter(paste(alpha, variable) != "0 lwr") %>% # remove duplicate in median (i.e. for alpha = 0, lwr = upr)
        mutate(quantile = paste0("Q", quantile*1000)) %>% 
        reshape2::dcast(vax_cov ~ quantile)
      if(plot_flag){plot_ind_boxcox_fits(
        error = model_errors, pred_nominal_pos = pred_nominal, to_plot = c("pos"), 
        plot_title = plot_title)}
    }
    else if(nrow(pos_errors) < fit_thresh){ # only fit negatives
      print("fitting only negatives")
      pred_nominal = fit_one_boxcox(
        neg_errors, new_vax_cov, alphas, neg_flag = TRUE
      )
      r <- pred_nominal %>%
        filter(paste(alpha, variable) != "0 lwr") %>% # remove duplicate in median (i.e. for alpha = 0, lwr = upr)
        mutate(quantile = paste0("Q", quantile*1000)) %>% 
        reshape2::dcast(vax_cov ~ quantile)
      if(plot_flag){plot_ind_boxcox_fits(
        error = model_errors, pred_nominal_neg = pred_nominal, to_plot = c("neg"), 
        plot_title = plot_title)}
    }
    else { # fit both
      print("fitting both")
      pred_nominal_pos = fit_one_boxcox(pos_errors, new_vax_cov, alphas) %>% 
        filter(paste(alpha, variable) != "0 lwr") # remove duplicate in median (i.e. for alpha = 0, lwr = upr)
      pred_nominal_neg = fit_one_boxcox(neg_errors, new_vax_cov, alphas, neg_flag = TRUE) %>%
        filter(paste(alpha, variable) != "0 lwr") # remove duplicate in median (i.e. for alpha = 0, lwr = upr)
      # combine into single distribution
      draws_neg = runif(ceiling(n_draws*prop_neg))
      draws_pos = runif(ceiling(n_draws*(1-prop_neg)))
      samp_summary = list()
      for(k in 1:nrow(new_vax_cov)){
        d_neg = pred_nominal_neg %>% 
          filter(vax_cov == new_vax_cov$vax_cov[k]) %>% 
          arrange(quantile)
        d_pos = pred_nominal_pos %>% 
          filter(vax_cov == new_vax_cov$vax_cov[k]) %>% 
          arrange(quantile)
        samp_neg = approx(x = d_neg$quant, y = d_neg$value, xout = draws_neg)$y
        samp_pos = approx(x = d_pos$quant, y = d_pos$value,xout = draws_pos)$y
        samp_summary[[k]] = data.frame(
          vax_cov = new_vax_cov$vax_cov[k], quantile = quantiles, 
          value = quantile(c(samp_neg, samp_pos), quantiles, na.rm = TRUE)
        )
      }
      r <- bind_rows(samp_summary) %>%
        mutate(quantile = paste0("Q", quantile*1000)) %>% 
        reshape2::dcast(vax_cov ~ quantile)
      if(plot_flag){plot_ind_boxcox_fits(
        error = model_errors, pred_nominal_pos = pred_nominal_pos, 
        pred_nominal_neg = pred_nominal_neg, to_plot = c("pos", "neg"), 
        plot_title = plot_title)}
    }
    return(r)
}

plot_ind_boxcox_fits <- function(error, pred_nominal_pos, pred_nominal_neg, 
                                 to_plot = c("neg", "pos"), plot_title){
  p <- ggplot(data = error, aes(x = vax_cov)) + 
    geom_point(aes(y = error)) + 
    ggtitle(plot_title) + 
    theme_bw()
  if("neg" %in% to_plot){
    p <- p + 
      geom_ribbon(data = pred_nominal_neg %>% filter(round(alpha,3) == 0.9) %>% 
                    reshape2::dcast(vax_cov + fit ~ variable, value.var = "value"),
                  aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "red") + 
      geom_ribbon(data = pred_nominal_neg %>% filter(round(alpha,3) == 0.5) %>% 
                    reshape2::dcast(vax_cov + fit ~ variable, value.var = "value"),
                  aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "red") + 
      geom_line(data = pred_nominal_neg %>% filter(round(quantile,3) == 0.9), 
                aes(y = fit), color = "red", size = 1)
  }
  if("pos" %in% to_plot){
    p <- p + 
      geom_ribbon(data = pred_nominal_pos %>% filter(round(alpha,3) == 0.9) %>% 
                    reshape2::dcast(vax_cov + fit ~ variable, value.var = "value"),
                  aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") + 
      geom_ribbon(data = pred_nominal_pos %>% filter(round(alpha,3) == 0.5) %>% 
                    reshape2::dcast(vax_cov + fit ~ variable, value.var = "value"),
                  aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") + 
      geom_line(data = pred_nominal_pos %>% filter(round(quantile,3) == 0.9), 
                aes(y = fit), color = "blue", size = 1) 
  }
  print(p)
}

run_full_analysis <- function(all_errors, new_vax_cov, alphas, n_draws = 1e5, fit_thresh = 10,
                              plot_ind_boxcox_fits_flag = FALSE, calculate_coverage_flag = TRUE, 
                              plot_full_results_flag = TRUE, full_results_title = NA){
  # subset to only observed errors (i.e., scenario_id = T returns only observed errors)
  observered_errors = all_errors %>% filter(scenario_id == "T")
  # subset to only true errors to test against (i.e., scenario_id == "S1", "S2")
  true_errors = all_errors %>% filter(scenario_id %in% c("S1", "S2"))
  # create container to store results
  full_results_list <- list()
  # fit mixture distribution for each mdoel
  preds_all <- vector("list", length = n_models)
  for(i in 1:n_models){
    print(i)
    dat_filt = observered_errors %>% filter(model_id == paste0("M", i))
    preds_all[[i]] = fit_mixture_w_boxcox(
      dat_filt, new_vax_cov, alphas, plot_flag = plot_ind_boxcox_fits_flag, 
      plot_title = paste("Model:", i)) %>% 
      mutate(model_id = paste0("M", i))
  }
  full_results_list$prediction = bind_rows(preds_all)
  # calculate coverage
  if(calculate_coverage_flag){
    cov = calculate_coverage(
      # do some reshaping to match expected format for quant_fits
      quant_fits = bind_rows(preds_all) %>% 
        melt(c("model_id", "vax_cov")) %>% 
        filter(variable != "mean") %>%
        mutate(quantile = as.integer(gsub("Q", "", variable))/1000) %>%
        dplyr::select(-variable), 
      error_df = true_errors,  # %>% mutate(error = abs(error))
      vax_cov_S1 = 0.3, vax_cov_S2 = 0.5, 
      summarize_by = "model"
    )
    full_results_list$coverage = cov
  }
  # plot results
  if(plot_full_results_flag){
    p1 = bind_rows(preds_all) %>% 
      ggplot(aes(x = vax_cov)) + 
      # geom_point(data = t_large$errors %>% filter(scenario_id %in% c("S1", "S2")),
      #            aes(y = error), color = "darkgray", shape = 8) +
      # geom_line(data = t_large$errors %>% filter(scenario_id == "E"),
      #           aes(y = error, group = location_id), color = "darkgray", alpha = 0.25) +
      geom_point(data = observered_errors, aes(y = error), shape = 21) +
      geom_ribbon(aes(ymin = Q25, ymax = Q975, fill = model_id), alpha = 0.4) +
      geom_ribbon(aes(ymin = Q250, ymax = Q750, fill = model_id), alpha = 0.6) +
      geom_line(aes(y = Q500, color = model_id), size = 1) + 
      facet_wrap(vars(model_id)) + 
      labs(x = "vaccination coverage", y = "absolute error", 
           title = full_results_title) + 
      theme_bw() + 
      theme(legend.position = "none")
    p2 = ggplot(data = cov, aes(x = alpha, y = cov, group = model_id)) + 
      geom_line(aes(color = model_id)) + 
      geom_abline(size = 1) + 
      facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs), ncol = 1) + 
      labs(x = "expected coverage", y = "actual coverage") + 
      theme_bw() + 
      theme(legend.position = "none")
    p = cowplot::plot_grid(p1, p2, rel_widths = c(0.7, 0.3), align = "v", axis = "t") 
    full_results_list$p = p
  }
  return(full_results_list)
}
  

#### TRY THE FULL ANALYSIS -----------------------------------------------------
r_large <- run_full_analysis(
  all_errors = t_large$errors, new_vax_cov = new_vax_cov, alphas = alphas, 
  plot_ind_boxcox_fits_flag = FALSE, calculate_coverage_flag = TRUE, 
  plot_full_results_flag = TRUE, 
  full_results_title = "500 locations"
)

r_large$p
ggsave("R/sim-experiment-final_size/mixture-distribution/500locations.pdf", width = 14, height = 6)

#### REPEAT WITH FEWER LOCATIONS -----------------------------------------------
t_small <- full_sim(n_locations = 50, n_models = n_models,
                    vax_cov_S1 = 0.3, vax_cov_S2 = 0.5,
                    R0_lwr = 2, R0_upr = 3,
                    model_bias_ind_sd = 0.05,
                    seed = seed_id, fit_outcomes = FALSE)

r_small <- run_full_analysis(
  all_errors = t_small$errors, new_vax_cov = new_vax_cov, alphas = alphas, 
  plot_ind_boxcox_fits_flag = FALSE, calculate_coverage_flag = TRUE, 
  plot_full_results_flag = TRUE, 
  full_results_title = "50 locations"
)

r_small$p
ggsave("R/sim-experiment-final_size/mixture-distribution/50locations.pdf", width = 14, height = 6)

#### REPEAT WITH LARGER CORRELATION IN OBSERVATIONS ----------------------------
t_corr <- full_sim(n_locations = 50, n_models = n_models,
                   vax_cov_S1 = 0.3, vax_cov_S2 = 0.5,
                   R0_lwr = 2, R0_upr = 3, cov_R0 = -0.25,
                   seed = seed_id, fit_outcomes = FALSE)

r_corr <- run_full_analysis(
  all_errors = t_corr$errors, new_vax_cov = new_vax_cov, alphas = alphas, 
  plot_ind_boxcox_fits_flag = FALSE, calculate_coverage_flag = TRUE, 
  plot_full_results_flag = TRUE, 
  full_results_title = "50 locations, high correlation between R0 and true vax (-0.5)"
)

r_corr$p
ggsave("R/sim-experiment-final_size/mixture-distribution/50locations_highcorr.pdf", width = 14, height = 6)

