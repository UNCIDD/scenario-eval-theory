library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)

source("R/simulation/0-parameters.R")
source("R/simulation/0-helper-functions.R")

method_comparison_results = readRDS("output/simulation/performance_evaluation_results_location_specific_sens.rds")

sim_out = readRDS("output/simulation/sim_out_sens.rds")

true_sims_df = bind_rows(lapply(sim_out, function(i){i$true_sims}), .id = "rep_id") %>%
  mutate(rep_id = as.integer(rep_id))

loc_labs =  true_sims_df %>%
  filter(scenario_id == "T") %>% 
  select(rep_id, location_id, location_R0, vax_cov) %>%
  rename(true_vax_cov = vax_cov) %>%
  unique()


ggplot(true_sims_df %>% filter(scenario_id == "T", rep_id %in% 1:10), 
       aes(x = vax_cov, y = true_final_size)) + 
  geom_point() + 
  
  facet_wrap(vars(rep_id))


obs_errors = approach3_obs_covariate %>% 
  filter(rep_id %in% 1:50) %>%
  reframe(est_obs_samp = get_samps(quantile, value, 1e2),
            draw_id = 1:1e2, .by = c("vax_cov", "scenario_id", "rep_id", "location_R0")) %>%
  left_join(true_sims_df) %>%
  mutate(obs_error = est_obs_samp - true_final_size) 
  

approach3_obs_covariate_ex = readRDS( "output/simulation/estimated_obs_approach3.rds")
sim_out_ex = readRDS("output/simulation/sim_out.rds")
obs_df_ex = sim_out_ex$true_sims

p1 = ggplot(obs_df_ex, aes(x = vax_cov, y = true_final_size, color = location_R0)) + 
  geom_point(data = obs_df_ex %>% filter(scenario_id == "T"), size = 2) + 
  # geom_line(aes(group = location_id)) +
  geom_line(data = approach3_obs_covariate_ex %>% filter(quantile == 0.5), 
            aes(x = vax_cov, y = value, group = location_id, color = location_R0), alpha = 0.5) +
  labs(x = "realized vaccine uptake\n(scenario axis)", y = "realized final size\n(observations)", color = "location\nR0") +
  scale_color_viridis_c() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.9)) +
  theme_bw() + 
  theme(panel.grid = element_blank())
p2 = ggplot(obs_df_ex %>% filter(location_id %in% c(44, 41))) + 
  # geom_line(aes(group = location_id)) +
  geom_line(data = approach3_obs_covariate_ex %>% filter(quantile == 0.5, location_id %in% c(44, 41)),
            aes(x = vax_cov, y = value, group = location_id, color = location_R0), linetype = "dotted") +
  geom_ribbon(data = approach3_obs_covariate_ex %>% 
                filter(quantile %in% c(0.05, 0.95), location_id %in% c(44, 41)) %>%
                mutate(quantile = paste0("Q", quantile*100)) %>% 
                dcast(vax_cov + location_id + location_R0 ~ quantile, value.var = "value"),
            aes(x = vax_cov, ymin = Q5, ymax = Q95, fill = location_R0), alpha = 0.4) +
  geom_line(aes(x = vax_cov, y = true_final_size), color = "black") + 
  facet_wrap(vars(paste0("R0 = ", round(location_R0, 2))), ncol = 1) +
  labs(x = "realized vaccine uptake\n(scenario axis)", y = "realized final size\n(projection axis)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.9)) +
  scale_color_viridis_c() + 
  scale_fill_viridis_c() + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.background = element_blank())

# bias in estimates of observations across replicates
p3 = obs_errors %>%
  mutate(R0_group = cut(location_R0, breaks = seq(1.975, 3.025, 0.05), labels = seq(2, 3, 0.05))) %>%
  ggplot(aes(x = as.factor(R0_group), y = obs_error))+ 
  geom_boxplot(aes(x = as.factor(R0_group), y = obs_error), alpha = 0.1) + 
  geom_hline(yintercept = 0, color = "red", lwd = 0.5) +
  scale_x_discrete(name = "location R0", breaks = seq(2, 3, 0.25)) +
  scale_y_continuous(name = "error in estimated\nobservations", limits = c(-1,1)*0.3) + 
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs)) + 
  theme_bw() +
  theme(strip.background = element_blank())

# approach2 vs. approach3 performance
p4 = method_comparison_results %>% 
  filter(approach != "1") %>%
  dcast(model_id + scenario_id + location_id + rep_id ~ approach, value.var = "mae") %>% 
  left_join(loc_labs) %>%
  mutate(R0_group = cut(location_R0, breaks = seq(1.975, 3.025, 0.05), labels = seq(2, 3, 0.05))) %>%
  ggplot() + 
  geom_boxplot(aes(x = as.factor(R0_group), y = `2-covariates` - `3-covariates`), alpha = 0.1) + 
  geom_label(data = data.frame(scenario_id = "S1", x = 3, y = c(-1, 1)*0.2, lab = c("approach 2\nmore accurate", "approach 3\nmore accurate")), 
            aes(x = x, y = y, label = lab),  color = "red", size = 2) +
  geom_hline(yintercept = 0, color = "red", lwd = 0.5) +
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs)) + 
  scale_color_viridis_c(name = "location R0") + 
  scale_x_discrete(name = "location R0", breaks = seq(2, 3, 0.25)) +
  scale_y_continuous(name = "performance difference \n(approach 2- approach 3)", limits = c(-1,1)*0.2) + 
  theme_bw() +
  theme(strip.background = element_blank())

plot_grid(
  plot_grid(p1, p2, rel_widths = c(0.65, 0.35), labels = c("A", "B")), 
  p3, 
  p4, 
  labels = c(NA, "C", "D"),
  ncol  = 1, 
  align = "h", axis = "lr", 
  rel_heights = c(0.4, 0.35, 0.35)
)
ggsave("output/figures/supp_location-specific-mae-a2vsa3_sens.pdf", width = 8, height = 10)
