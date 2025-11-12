library(ggplot2)
library(dplyr)

#### SETUP ---------------------------------------------------------------------
sim_out = readRDS("output/simulation/sim_out.rds")
all_ests_all_locs = readRDS("output/simulation/estimated_errors_across_locations.rds")

approach2_nocov_errors_across_locs = readRDS("output/simulation/estimated_errors_approach2_full_across_locs.rds")

error_df_small = sim_out$errors %>% filter(scenario_id == "T")
approach2_cov_ests = approach2_nocov_errors_across_locs %>% 
  filter(round(quantile,4) %in% c(0.05, 0.25, 0.5, 0.75, 0.95)) %>%
  mutate(quantile = paste0("Q", quantile*100)) %>% 
  dcast(vax_cov + scenario_id + model_id ~ quantile)

### APPROACH 2 FIGURE ----------------------------------------------------------
locs_to_plot = c(29, 5, 45)
loc_labs = c(15, 41, 32)
loc_cols = RColorBrewer::brewer.pal(4, "Set1")[2:4]
model_to_plot = "M8"

# overall results
mod_results = ggplot(data = approach2_cov_ests %>% filter(model_id == model_to_plot), 
                     aes(x = vax_cov)) +
  geom_point(data = error_df_small %>% filter(model_id == model_to_plot), 
             aes(x = vax_cov, y = error), color = "darkgray", size = 2, shape = 21) +
  geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.15, fill = "darkgray") +
  geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.3, fill = "darkgray") +
  geom_line(aes(y = Q50), alpha = 0.8, size = 1.5, color = "darkgray") +
  geom_point(data = error_df_small %>% filter(model_id == model_to_plot, location_id %in% locs_to_plot), 
             aes(y = error, color = as.factor(location_id)), size = 2) + 
  geom_text(data = error_df_small%>% 
              filter(model_id == model_to_plot, location_id %in% locs_to_plot) %>% 
              mutate(loc_lab = loc_labs[which(locs_to_plot == location_id)], .by = "location_id"),
            aes(y = error, label = loc_lab), size = 1.4, color = "white") +
  facet_wrap(vars(model_id), scales = "free") +
  labs(x = "realized vaccine uptake\n(scenario axis)", 
       subtitle = "Step 2: fit errors across scenario axis and infer error in modeled scenarios") +
  scale_color_manual(values = loc_cols) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", 
        panel.grid.minor =  element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank())
mod_results

# a few specific locations
loc_results = ggplot(data = sim_out$model_sims %>% 
                       filter(model_id %in% c(model_to_plot, "T"), location_id %in% locs_to_plot, 
                              scenario_id != "E") %>%
                       mutate(location_id = factor(location_id, levels = locs_to_plot))) + 
  geom_line(data = sim_out$model_sims %>% 
              filter(model_id == model_to_plot, location_id %in% locs_to_plot) %>%
              mutate(location_id = factor(location_id, levels = locs_to_plot)), 
            aes(x = vax_cov, y = final_size), alpha = 0.2) + 
  geom_point(aes(x = vax_cov, y = final_size, shape = scenario_id, fill = scenario_id), color = 'black', size = 1.5) +
  geom_text(data = data.frame(location_id = locs_to_plot, 
                              location_lab = loc_labs) %>%
              mutate(location_id = factor(location_id, levels = locs_to_plot)),
            aes(x = Inf, y = Inf, label = paste0("location ", location_lab), color = as.factor(location_id)), 
            hjust = 1, vjust = 1, size = 2.5) +
  geom_point(data = sim_out$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id == "T"), 
             aes(x = vax_cov, y = true_final_size), color = "red", size = 1.5) + 
  geom_segment(data = left_join(
    sim_out$model_sims %>% filter(model_id == model_to_plot, location_id %in% locs_to_plot, scenario_id == "T"),
    sim_out$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id == "T")
  ), aes(x = vax_cov, xend = vax_cov, y = true_final_size, yend = final_size), arrow = arrow(length = unit(0.05, "npc")), linewidth = 0.3, size = 0.4) +
  facet_wrap(vars(location_id), ncol = 1) +
  labs(x = "realized vaccine uptake\n(scenario axis)", 
       y = "final epidemic size\n(projection axis)", 
       subtitle = "Step 1: calculate error in realized scenario") +
  scale_color_manual(values = loc_cols[sapply(sort(locs_to_plot), function(i){which(locs_to_plot == i)})]) +
  scale_fill_manual(values = c("black", "black", "white")) +
  scale_shape_manual(values = c(16, 16, 21)) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank())
loc_results

cowplot::plot_grid(loc_results, mod_results, rel_widths = c(0.35, 0.65), 
                   align = "h", axis = "tb",
                   labels = c("A", "B"), label_size = 10)

ggsave("output/figures/approach2_illustration.pdf", width = 6.45, height = 3.5)

