library(ggplot2)
library(dplyr)
library(cowplot)

#### SETUP ---------------------------------------------------------------------
source("R/simulation/0-parameters.R")
source("R/simulation/0-helper-functions.R")

sim_out = readRDS("output/simulation/sim_out.rds")
all_ests = readRDS("output/simulation/estimated_errors_across_locations.rds")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#### GET INFERRED DISTRIBUTION OF OBSERVATIONS ---------------------------------
inferred_obs = all_ests %>%
  filter(substr(approach,1,1) == "2") %>% 
  # get samples of miscalibration error
  reframe(miscal_error = get_samps(quantile, value, 1e4),
                  draw_id = 1:1e4, .by = c("vax_cov", "scenario_id", "model_id", "approach")
  ) %>%
  # get projected value
  left_join(sim_out$model_sims %>% filter(scenario_id %in% c("S1", "S2")) %>% select(model_id, location_id, scenario_id, vax_cov, final_size), 
            relationship = "many-to-many") %>%
  # calculate inferred observation
  mutate(inferred_obs = final_size - miscal_error)  %>% 
  mutate(model_id = factor(model_id, levels = paste0("M", 1:10)))
# do KS test on these dists too
ks_test_rslts = expand.grid(model_id = paste0("M", 1:n_models), 
                            scenario_id = c("S1", "S2"))
ks_test_rslts$ks_test_stat = NA
ks_test_rslts$ks_test_p = NA
ks_test_rslts$n_df = NA
ks_test_rslts$n_truth = NA
for(i in 1:nrow(ks_test_rslts)){
  tmp_scenario_id = ks_test_rslts[i, "scenario_id"]
  tmp_model_id = ks_test_rslts[i, "model_id"]
  # pull error distribution of interest
  tmp_df = inferred_obs %>% filter(scenario_id == tmp_scenario_id, model_id == tmp_model_id)
  # get true error
  tmp_truth = sim_out$true_sims %>% filter(scenario_id == tmp_scenario_id)
  # perform KS test
  tmp_ks = ks.test(tmp_df$inferred_obs, tmp_truth$true_final_size)
  # save output
  ks_test_rslts[i, "n_df"] = length(tmp_df$inferred_obs)
  ks_test_rslts[i, "n_truth"] = length(tmp_truth$true_final_size)
  ks_test_rslts[i, "ks_test_stat"] = tmp_ks$statistic
  ks_test_rslts[i, "ks_test_p"] = tmp_ks$p.value
}

p1 = ggplot(data = inferred_obs, 
       aes(x = inferred_obs, color = model_id)) + 
  geom_density() + 
  geom_density(data = sim_out$true_sims %>% filter(scenario_id %in% c("S1", "S2")), 
               aes(x = true_final_size, color = "actual")) + 
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs), ncol = 1) +
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual(values = c(gg_color_hue(10), "black")) + 
  scale_x_continuous(name = "observed final size") + 
  theme_bw() + 
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())

p2 = ks_test_rslts %>% 
  ggplot(aes(x = model_id, y = ks_test_stat, shape = scenario_id, color = model_id)) + 
  geom_hline(data = ks_test_rslts %>%
               dplyr::select(n_df, n_truth) %>%
               unique() %>%
               mutate(ks_sig = 1.358*sqrt((n_df + n_truth)/(n_df*n_truth))), # 5% p value level 
             aes(yintercept = ks_sig), linewidth = 0.3) + 
  geom_point(size = 2.5) + 
  guides(color = "none") +
  scale_shape_discrete(labels = scenario_labs) + 
  scale_y_continuous(name = "Komogorov-Smirnov test statistic") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        panel.grid.major.x = element_blank())

plot_grid(p1, p2, labels = c("A", "B"), align = "h", axis = "tb",
          nrow = 1, rel_widths = c(0.4, 0.6))

ggsave("output/figures/supp_inferred-error-dists-across-locations.pdf", width = 8, height = 5)

#### REPEAT WITH COVARIATES ----------------------------------------------------
all_ests = readRDS("output/simulation/estimated_errors_location_specific.rds")

loc_labs = sim_out$true_sims %>%
  filter(scenario_id == "T") %>% 
  select(location_id, location_R0) %>%
  unique() %>%
  mutate(lab = paste0("location ", location_id, " (R0 = ", round(location_R0, 2), ")"), 
         lab2 = paste0("R0 = ", round(location_R0, 2), " (location ", location_id, ")"))
loc_labs = loc_labs %>% 
  mutate(lab = factor(lab, levels = loc_labs$lab), 
         lab2 = factor(lab2, levels = sort(loc_labs$lab2)))


inferred_obs = all_ests %>%
  filter(substr(approach,1,1) == "2") %>% 
  # get samples of miscalibration error
  reframe(miscal_error = get_samps(quantile, value, 1e4),
          draw_id = 1:1e4, .by = c("vax_cov", "scenario_id", "location_id", "model_id", "approach")
  ) %>%
  # get projected value
  left_join(sim_out$model_sims %>% filter(scenario_id %in% c("S1", "S2")) %>% select(model_id, location_id, scenario_id, vax_cov, final_size)) %>%
  # calculate inferred observation
  mutate(inferred_obs = final_size - miscal_error)  %>%
  mutate(model_id = factor(model_id, levels = paste0("M", 1:10)))

ggplot(data = inferred_obs %>% left_join(loc_labs), 
       aes(x = inferred_obs, color = model_id, group = interaction(model_id, scenario_id))) + 
  geom_density(linewidth = 0.3) + 
  geom_vline(data = sim_out$true_sims %>% filter(scenario_id %in% c("S1", "S2")) %>% left_join(loc_labs), 
             aes(xintercept = true_final_size, linetype = scenario_id, color = "actual")) + 
  facet_wrap(vars(lab2)) +
  scale_color_manual(values = c(gg_color_hue(10), "black")) + 
  scale_linetype_manual(values = c("solid", "longdash"), labels = scenario_labs) + 
  scale_x_continuous(name = "final size in modeled scenarios", limits = c(0, 1)) + 
  theme_bw(base_size = 7) + 
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank())
ggsave("output/figures/supp_inferred-error-dists-location-specific.pdf", width = 8, height = 5)

