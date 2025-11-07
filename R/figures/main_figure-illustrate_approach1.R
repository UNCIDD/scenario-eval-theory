library(ggplot2)
library(dplyr)


#### SETUP ---------------------------------------------------------------------
sim_out = readRDS("output/simulation/sim_out.rds")

### APPROACH 1 FIGURE ----------------------------------------------------------
locs_to_plot = c(29, 13, 5) #19, 
loc_labs = c(15, 39, 41)
loc_plaus = c(rep("S1", 2), rep("S2", 1))
loc_cols = RColorBrewer::brewer.pal(4, "Set1")[2:4]
model_to_plot = "M8"

approach_1_plot_df = sim_out$model_sims %>% 
  filter(model_id %in% c(model_to_plot, "T"), location_id %in% locs_to_plot, 
         scenario_id %in% c("S1", "S2")) %>%
  left_join(sim_out$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id == "T") %>% rename(true_vax_cov = vax_cov) %>% dplyr::select(-scenario_id)) %>%
  left_join(data.frame(location_id = locs_to_plot, plausible_scenario = loc_plaus)) %>%
  mutate(location_id = factor(location_id, levels = locs_to_plot))

ggplot(data = approach_1_plot_df) + 
  geom_line(data = sim_out$true_sims %>% filter(location_id %in% locs_to_plot), 
            aes(x = vax_cov, y = true_final_size, color = "observation"), alpha = 0.2, linewidth = 0.5) +
  geom_segment(data = sim_out$model_sims %>% 
                 filter(model_id %in% c(model_to_plot, "T"), location_id %in% locs_to_plot, 
                        scenario_id %in% c("S1", "S2")) %>%
                 left_join(sim_out$true_sims %>% 
                             filter(location_id %in% locs_to_plot, 
                                    scenario_id %in% c("S1", "S2"))) %>%
                 mutate(location_id = factor(location_id, levels = locs_to_plot)), 
               aes(x = vax_cov,# + ifelse(scenario_id == "S1", -0.005, 0.005), 
                   xend = vax_cov,# + ifelse(scenario_id == "S1", -0.005, 0.005), 
                   y = final_size, yend = true_final_size), linewidth = 1, alpha = 0.2) + 
  geom_point(data = sim_out$true_sims %>% filter(location_id %in% locs_to_plot, scenario_id %in% c("S1", "S2")), 
             aes(x = vax_cov, y = true_final_size, color = "observation"), shape = 21, fill = "white", alpha = 0.2) +
  geom_segment(data = approach_1_plot_df %>% filter(scenario_id == plausible_scenario),
               aes(x = vax_cov, xend = true_vax_cov, y = true_final_size, yend = true_final_size), linetype = "dotted", linewidth = 0.3) + 
  geom_segment(data = approach_1_plot_df %>% filter(scenario_id == plausible_scenario),
               aes(x = vax_cov, xend = vax_cov, 
                   y = true_final_size, yend = final_size), arrow = arrow(length = unit(0.035, "npc")), linewidth = 0.3) +
  geom_point(aes(x = vax_cov, y = final_size, color = "projection"), size = 1.25) + 
  geom_point(aes(x = true_vax_cov, y = true_final_size, color = "observation"), size = 1.25) + 
  geom_text(data = data.frame(location_id = locs_to_plot, 
                              location_lab = loc_labs) %>%
              mutate(location_id = factor(location_id, levels = locs_to_plot)),
            aes(x = Inf, y = Inf, label = paste0("location ", location_lab)), 
            hjust = 1, vjust = 1, size = 2.5) +
  facet_wrap(vars(location_id)) +
  labs(x = "realized vaccine uptake\n(scenario axis)", 
       y = "cumulative hospitalizations\n(projection axis)") +
  scale_color_manual(values = c("red", "black")) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank())
ggsave("output/figures/approach1_illustration.pdf", width = 5, height = 2.5)

