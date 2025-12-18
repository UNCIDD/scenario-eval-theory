library(ggplot2)
library(dplyr)

#### SETUP ---------------------------------------------------------------------
sim_out = readRDS("output/simulation/sim_out.rds")
approach3_obs_covariate = readRDS("output/simulation/estimated_obs_approach3.rds")

approach3_cov_ests = approach3_obs_covariate %>%
  filter(round(quantile,4) %in% c(0.05, 0.25, 0.5, 0.75, 0.95)) %>%
  mutate(quantile = paste0("Q", quantile*100)) %>% 
  dcast(vax_cov + scenario_id + location_id ~ quantile)

### APPROACH 3 FIGURE ----------------------------------------------------------
locs_to_plot = c(31, 45) #29
loc_labs = c(9, 32) #15
loc_cols = rev(RColorBrewer::brewer.pal(5, "Set1")[4:5])
model_to_plot = "M8"

plot_fit_obs = ggplot(data = sim_out$true_sims %>% filter(scenario_id == "T"), 
                      aes(x = vax_cov)) + 
  geom_ribbon(data = approach3_cov_ests %>% filter(location_id %in% locs_to_plot), 
              aes(ymin = Q5, ymax = Q95, fill = as.factor(location_id)), alpha = 0.25) +
  # geom_line(data = bind_rows(fit_obs_noR0) %>% filter(alpha == 0.95) %>% 
  #               mutate(lwr = ifelse(is.na(lwr), 0, lwr)), 
  #             aes(y = fit), alpha = 0.8, color = "red", linetype = "dashed") +
  geom_line(data = approach3_cov_ests %>% filter(location_id %in% locs_to_plot), 
            aes(y = Q50, color = as.factor(location_id)), linewidth = 0.3, linetype = "dashed") +
  geom_point(data = approach3_cov_ests %>% filter(location_id %in% locs_to_plot, vax_cov %in% c(0.3, 0.5)), 
             aes(y = Q50, color = as.factor(location_id)),shape = 21, fill = "white", size = 1.5) +
  geom_point(aes(y = true_final_size), color = "red", size = 1.5) + 
  geom_point(data = sim_out$true_sims %>% filter(scenario_id == "T", location_id %in% locs_to_plot), 
             aes(y = true_final_size, color = as.factor(location_id)), size = 1.5) +
  geom_text(data = sim_out$true_sims %>% filter(scenario_id == "T") %>%
              filter(location_id %in% locs_to_plot) %>% 
              left_join(data.frame(location_id = locs_to_plot, 
                                   new_location_id = loc_labs)), 
            aes(y = true_final_size, label = new_location_id), color = "white", size = 1) + 
  labs(x = "realized vaccine uptake\n(scenario axis)", 
       y = "observed final epidemic size\n(projection axis)", 
       subtitle = "Step 1: fit observations across scenario axis") +
  scale_color_manual(values = loc_cols) +
  scale_fill_manual(values = loc_cols) + 
  # scale_y_continuous(limits = c(0.1, 0.8)) + 
  theme_bw(base_size = 7) + 
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())
plot_fit_obs

approach_3_plot_df = sim_out$model_sims %>% 
  filter(model_id %in% c(model_to_plot, "T"), location_id %in% locs_to_plot, 
         scenario_id %in% c("S1", "S2")) %>%
  left_join(sim_out$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id %in% c("S1", "S2"))) %>%
  mutate(location_id = factor(location_id, levels = locs_to_plot))


plot_calc_error = ggplot(data = approach_3_plot_df, aes(x = vax_cov)) + 
  geom_line(data = sim_out$true_sims %>% filter(location_id %in% locs_to_plot), 
            aes(y = true_final_size), linewidth = 0.5, color = "red", alpha = 0.3) + 
  geom_ribbon(data = approach3_cov_ests %>% filter(location_id %in% locs_to_plot), 
              aes(ymin = Q5, ymax = Q95, fill = as.factor(location_id)), alpha = 0.25) +
  geom_line(data = approach3_cov_ests %>% filter(location_id %in% locs_to_plot), 
            aes(y = Q50, color = as.factor(location_id)), linewidth = 0.3, linetype = "dashed") +
  geom_point(data = approach3_cov_ests %>% filter(location_id %in% locs_to_plot, vax_cov %in% c(0.3, 0.5)), 
             aes(y = Q50, color = as.factor(location_id)), shape = 21, fill = "white", size = 1.5) +
  geom_point(aes(y = final_size), size = 1.5) + 
  geom_point(aes(y = true_final_size), size = 1.5, shape = 21, color = "red", fill = "white", alpha = 0.3) +
  geom_segment(data = approach3_cov_ests %>%
                 filter(location_id %in% locs_to_plot, vax_cov %in% c(0.3, 0.5)) %>%
                 left_join(sim_out$model_sims %>%
                             filter(model_id %in% c(model_to_plot, "T"), location_id %in% locs_to_plot,
                                    scenario_id %in% c("S1", "S2"))),
               aes(x = vax_cov, xend = vax_cov,
                   y = Q50, yend = final_size), arrow = arrow(length = unit(0.035, "npc")), linewidth = 0.3, size = 0.4) +
  geom_text(data = data.frame(location_id = locs_to_plot, 
                              location_lab = loc_labs) %>%
              mutate(location_id = factor(location_id, levels = locs_to_plot)),
            aes(x = Inf, y = Inf, label = paste0("location ", location_lab), color = as.factor(location_id)), 
            hjust = 1, vjust = 1, size = 2.5) +
  facet_wrap(vars(location_id), ncol = 1, scales = "free") +
  scale_color_manual(values = loc_cols) +
  scale_fill_manual(values = loc_cols) + 
  # scale_y_continuous(limits = c(0.1, 0.8)) + 
  labs(x = "realized vaccine uptake\n(scenario axis)", 
       y = "final epidemic size\n(projection axis)", 
       subtitle = "Step 2: calculate error from inferred observations") +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", 
        legend.title = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank())
plot_calc_error

cowplot::plot_grid(plot_fit_obs, plot_calc_error, 
                   rel_widths = c(0.575, 0.425), 
                   nrow = 1, labels = c("A", "B"), label_size = 10)  
ggsave("output/figures/approach3_illustration.pdf", width = 6.25, height = 3.5)

