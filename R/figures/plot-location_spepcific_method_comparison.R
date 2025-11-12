library(ggplot2)
library(dplyr)
library(cowplot)

#### SETUP ---------------------------------------------------------------------
source("R/simulation/0-parameters.R")

sim_out = readRDS("output/simulation/sim_out.rds")
all_ests = readRDS("output/simulation/estimated_errors_location_specific.rds")
method_comparison_results = readRDS("output/simulation/performance_evaluation_results_location_specific.rds")
  

loc_labs = sim_out$true_sims %>%
  filter(scenario_id == "T") %>% 
  select(location_id, location_R0, vax_cov) %>%
  rename(true_vax_cov = vax_cov) %>%
  unique() %>%
  mutate(lab = paste0("location ", location_id, " (R0 = ", round(location_R0, 2), ")"), 
         lab2 = paste0("R0 = ", round(location_R0, 2), " (location ", location_id, ")"))
loc_labs = loc_labs %>% 
  mutate(lab = factor(lab, levels = loc_labs$lab))

#### PLOT RESULTS --------------------------------------------------------------

# mae for all locations
p1 = ggplot(data = method_comparison_results) + 
  geom_point(aes(x = substr(approach, 1, 1), y = mae, color = approach), alpha = 0.2, position = "dodge") + 
  geom_violin(aes(x = substr(approach, 1, 1), y = mae, color = approach), draw_quantiles = 0.5, fill = NA) +
  facet_grid(rows = vars(scenario_id), cols = vars(model_id), 
             labeller = labeller(scenario_id = scenario_labs), switch = "y") + 
  scale_color_brewer(palette = "Set1") + 
  scale_x_discrete(name = "estimation approach") + 
  scale_y_continuous(name = "absolute difference between means") +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit(0, "lines"),
        strip.background = element_blank(), 
        strip.placement = "outside")

p2 = ggplot(data = method_comparison_results) + 
  geom_point(aes(x = substr(approach, 1, 1), y = mae, color = approach), alpha = 0.5, position = "dodge") + 
  geom_violin(aes(x = substr(approach, 1, 1), y = mae, color = approach), draw_quantiles = 0.5, fill = NA) +
  facet_grid(cols = vars(scenario_id), labeller = labeller(scenario_id = scenario_labs)) + 
  scale_color_brewer(palette = "Set1") + 
  scale_x_discrete(name = "estimation approach") + 
  scale_y_continuous(name = "absolute difference between means") +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_blank())
plot_grid(p1, p2, ncol = 1, labels = c("A", "B"))
ggsave("output/figures/supp_location-specific-mae.pdf", width = 12, height = 8)

ggplot(data = method_comparison_results %>% 
         filter(approach != "1") %>%
         dcast(model_id + scenario_id + location_id ~ approach, value.var = "mae") %>% 
         left_join(loc_labs)) + 
  geom_abline() + 
  geom_point(aes(x = `2-covariates`, y = `3-covariates`, color = location_R0), size = 2, alpha = 0.5) + 
  facet_wrap(vars(scenario_id), labeller = labeller(scenario_id = scenario_labs)) + 
  scale_color_viridis_c(name = "location R0") + 
  scale_x_continuous(name = "approach 2\nabsolute difference between means") + 
  scale_y_continuous(name = "approach 3\nabsolute difference between means") + 
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("output/figures/supp_location-specific-mae-a2vsa3.pdf", width = 8, height = 4)

method_comparison_results %>% 
  left_join(loc_labs) %>% 
  mutate(scenario_vax_cov = ifelse(scenario_id == "S1", vax_scenarios[1], vax_scenarios[2])) %>%
  ggplot(aes(x = scenario_vax_cov - true_vax_cov, y = mae, color = scenario_id)) + 
  geom_point() +
  facet_wrap(vars(approach), labeller = labeller(approach = approach_labs)) + 
  scale_color_manual(values = c("black", "gray"), labels = scenario_labs) +
  scale_x_continuous(name = "(scenario assumed coverage) - (realized coverage)") + 
  scale_y_continuous(name = "absolute difference between means") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        strip.background = element_blank())
ggsave("output/figures/supp_location-specific-mae-byscenarval.pdf", width = 9, height = 4)

