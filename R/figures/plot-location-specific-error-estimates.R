#### SETUP ---------------------------------------------------------------------
all_ests_all_locs = readRDS("output/simulation/estimated_errors_location_specific.rds")

# setup plotting
set.seed(68)
locations_to_plot = sample(1:n_loc, 6)
models_to_plot = paste0("M", sample(1:n_models, 2))

loc_labs = obs_df %>% 
  select(location_id, location_R0) %>%
  unique() %>%
  mutate(lab = paste0("location ", location_id, " (R0 = ", round(location_R0, 2), ")"))
loc_labs = loc_labs %>% 
  mutate(lab = factor(lab, levels = loc_labs$lab))

# plot estimated errors for each location (when possible) 
all_ests_all_locs %>%
  filter(round(quantile,2) %in% c(0.05, 0.25, 0.5, 0.75, 0.95), scenario_id %in% c("S1", "S2"), 
         location_id %in% locations_to_plot, model_id %in% models_to_plot, approach != "2-nocovariates") %>%
  mutate(quantile = paste0("Q", quantile*100)) %>%
  dcast(vax_cov + scenario_id + model_id + location_id + approach ~  quantile, value.var = "value") %>%
  left_join(loc_labs) %>%
  mutate(approach = factor(approach, levels = c("truth", "1", "2-covariates", "3-covariates")), 
         scenario_id = factor(scenario_id, levels = c("S2", "S1"))) %>%
  mutate(yval = as.integer(as.factor(scenario_id)) + as.integer(as.factor(approach))*-0.2 - 0.5) %>%
  ggplot() + 
  geom_hline(yintercept = 0.5, color = "gray", size = 0.3) + 
  geom_segment(data = all_ests_all_locs %>%
                 filter(quantile == 0.5, scenario_id %in% c("S1", "S2"), 
                        location_id %in% locations_to_plot, model_id %in% models_to_plot, approach == "truth") %>%
                 mutate(scenario_id = factor(scenario_id, levels = c("S2", "S1"))) %>% 
                 mutate(ystart = as.integer(as.factor(scenario_id)) - 1.4, yend = as.integer(as.factor(scenario_id)) - 0.7) %>%
                 left_join(loc_labs), 
               aes(x = value, xend = value, y = ystart, yend = yend), alpha = 0.5, size = 0.5) + 
  geom_segment(aes(x = Q5, xend = Q95, y = yval, yend = yval, color = approach), size = 0.6) + 
  geom_segment(aes(x = Q25, xend = Q75, y = yval, yend = yval, color = approach), size = 1) + 
  geom_point(aes(x = Q50, y = yval, color = approach)) + 
  facet_grid(cols = vars(lab), rows = vars(model_id), switch = "y", 
             labeller = labeller(model_id = model_labs)) + 
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(6, "Paired")[c(2, 4, 6)]), 
                     labels = approach_labs) + 
  scale_x_continuous("miscalibration error") + 
  scale_y_continuous(breaks = 0:1, labels = c("scenario 2\n(high vax)", "scenario 1\n(low vax)")) +
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        strip.placement = "outside")
ggsave("output/figures/supp_location-specific-error-ests.pdf", width = 12, height = 6)

