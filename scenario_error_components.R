library(ggplot2)
library(cowplot)

# set up simple example assuming linearity
ex_scenarios = c(0, 1)
obs_scenario = 0.7
ex_projections = c(140, 60)
intermed_projection = max(ex_projections) + diff(ex_projections)*obs_scenario
true_projections = c(90, 40)
true_obs = max(true_projections) + diff(true_projections)*obs_scenario

example_data = data.frame(
  scenarios = rep(c(ex_scenarios, obs_scenario),2), 
  projections = c(ex_projections, intermed_projection, true_projections, true_obs),
  type = c("model", "model", "model", "observation", "observation", "observation"), 
  observable = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
)

p_base <- ggplot(data = example_data) + 
  geom_vline(xintercept = ex_scenarios) + 
  geom_vline(xintercept = obs_scenario, linetype = "dotted") + 
  geom_line(aes(x = scenarios, y = projections, color = type), alpha = 0.3) + 
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 3) +
  guides(fill = FALSE) +
  scale_color_manual(values = c("black", "#ff420e"), 
                     labels = c("projection", "observation")) + 
  scale_fill_manual(values = c("white", "white", "black", "#ff420e")) + 
  scale_x_continuous(breaks = c(ex_scenarios, obs_scenario), 
                     labels = c("scenario 1", "scenario 2", "realized"), 
                     limits = c(- 0.1, 1.1),
                     name = "scenario axis") + 
  scale_y_continuous(name = "projection axis") + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid = element_blank())
  
p_miscal = p_base + 
  geom_segment(data = data.frame(x = obs_scenario, xend = obs_scenario, 
                                 y = intermed_projection, yend = true_obs),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#004586") + 
  geom_text(data = data.frame(x = obs_scenario-0.02, 
                              y = mean(c(intermed_projection, true_obs))),
            aes(x = x, y = y), label = "miscalibration\nbias", 
            color = "#004586", hjust = 1, size = 2) + 
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 3)

p_tot = p_base +   
  geom_segment(data = data.frame(x = ex_scenarios[1], xend = ex_scenarios[2], 
                                 y = true_obs, yend = true_obs),
               aes(x = x, xend = xend, y = y, yend = yend),  linetype = "dotted") + 
  geom_segment(data = data.frame(x = ex_scenarios[1], xend = ex_scenarios[1], 
                                 y = ex_projections[1], yend = true_obs),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#579d1c") + 
  geom_text(data = data.frame(x = ex_scenarios[1]+0.03, 
                              y = mean(c(ex_projections[1], true_obs))), 
            aes(x = x, y = y), label = "observation\nbias", 
            color = "#579d1c", hjust = 0, size = 2)  +
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 3)

p_true = p_base +   
  geom_segment(data = data.frame(x = ex_scenarios[1], xend = ex_scenarios[1],
                                 y = ex_projections[1], yend = true_projections[1]),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#83caff") +
  geom_segment(data = data.frame(x = ex_scenarios[2], xend = ex_scenarios[2],
                                 y = ex_projections[2], yend = true_projections[2]),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#83caff") +
  geom_text(data = data.frame(x = ex_scenarios[1]+0.03, 
                              y = mean(c(ex_projections[1], true_projections[1]))),
            aes(x = x, y = y),
            label = "true\nbias", color = "#83caff", hjust = 0, size = 2) + 
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 3)

 
p_scenario = p_base +   
  geom_segment(data = data.frame(x = ex_scenarios[1], xend = ex_scenarios[2], 
                                 y = intermed_projection, yend = intermed_projection),
               aes(x = x, xend = xend, y = y, yend = yend), linetype = "dotted") + 
  geom_segment(data = data.frame(x = ex_scenarios[1], xend = ex_scenarios[1],
                                 y = ex_projections[1], yend = intermed_projection),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#7e0021") +
  geom_segment(data = data.frame(x = ex_scenarios[2], xend = ex_scenarios[2],
                                 y = ex_projections[2], yend = intermed_projection),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#7e0021") +
  geom_text(data = data.frame(x = ex_scenarios[1]+0.03,
                              y = mean(c(ex_projections[1], intermed_projection))),
            aes(x = x, y = y), label = "scenario\nbias",
            color = "#7e0021", hjust = 0, size = 2) +
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 3)
 

plot_grid(p_true, p_tot, p_miscal, p_scenario)
ggsave("scenario_error_components.pdf", width = 6.5, height = 6)


