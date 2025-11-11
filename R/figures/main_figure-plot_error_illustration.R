library(ggplot2)
library(cowplot)
library(dplyr)

# set up simple example assuming linearity
ex_scenarios = c(0, 1)
obs_scenario = 0.7
# ex_projections = c(140, 0)
ex_projections = 140*0.38^ex_scenarios
intermed_projection = 140*0.38^obs_scenario
extra_proj = 140*(0.38^seq(ex_scenarios[1], ex_scenarios[2], length.out = 25))
true_projections = c(90, 40)
b = ((true_projections[2]/true_projections[1])^(1 /(ex_scenarios[2] - ex_scenarios[1])))
true_obs = true_projections[1]*(b^obs_scenario)
# true_obs = max(true_projections) + diff(true_projections)*obs_scenario
extra_obs = true_projections[1]*(b^seq(ex_scenarios[1], ex_scenarios[2], length.out = 25))


example_data = data.frame(
  scenarios = rep(c(ex_scenarios, obs_scenario),2), 
  projections = c(ex_projections, intermed_projection, true_projections, true_obs),
  type = c("model", "model", "model", "observation", "observation", "observation"), 
  observable = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
)

extra_data = bind_rows(
  data.frame(
    scenarios = seq(ex_scenarios[1], ex_scenarios[2], length.out = 25), 
    projections = extra_obs, 
    type = "observation"),
  data.frame(
    scenarios = seq(ex_scenarios[1], ex_scenarios[2], length.out = 25), 
    projections = extra_proj, 
    type = "model"),
)

p_base <- ggplot(data = example_data) + 
  geom_vline(xintercept = ex_scenarios, linewidth = 0.3) + 
  geom_vline(xintercept = obs_scenario, linetype = "dotted", linewidth = 0.3) + 
  geom_line(data = extra_data, aes(x = scenarios, y = projections, color = type), alpha = 0.3, linewidth = 0.5) +
  geom_segment(data = example_data %>% filter(scenarios != obs_scenario) %>% dcast(scenarios ~ type, value.var = "projections"), 
               aes(x = scenarios, xend = scenarios, y = model, yend = observation), linewidth = 1.25, color = "#83caff") + 
  geom_segment(data = example_data %>% filter(scenarios == obs_scenario) %>% dcast(scenarios ~ type, value.var = "projections"), 
               aes(x = scenarios, xend = scenarios, y = model, yend = observation), linewidth = 1.25, color = "#579d1c") + 
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 2) +
  guides(fill = "none") +
  scale_color_manual(values = c("black", "#ff420e"), 
                     labels = c("projection", "observation")) + 
  scale_fill_manual(values = c("white", "white", "black", "#ff420e")) + 
  scale_x_continuous(breaks = c(ex_scenarios, obs_scenario), 
                     labels = c("low vax", "high vax", "realized"), 
                     limits = c(- 0.1, 1.1),
                     name = "scenario axis") + 
  scale_y_continuous(limits = range(example_data$projections), 
                     name = "cumulative hospitalizations\n(projection axis)") + 
  theme_classic(base_size = 7) +
  theme(legend.position = "none", 
        panel.grid = element_blank())
plot_grid(p_base, labels = "A", label_size = 8)
ggsave("output/figures/base_plot.pdf", width = 2.5, height = 2.5)

