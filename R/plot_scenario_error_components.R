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
  geom_vline(xintercept = ex_scenarios) + 
  geom_vline(xintercept = obs_scenario, linetype = "dotted") + 
  geom_line(data = extra_data, aes(x = scenarios, y = projections, color = type), alpha = 0.3) +
  geom_point(aes(x = scenarios, y = projections, color = type, 
                 fill = interaction(type,observable)), shape = 21, size = 3) +
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
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid = element_blank())
p_base

ggsave("base_plot.tiff", p_base, width = 3.25, height = 3)
  
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
  geom_segment(data = data.frame(x = ex_scenarios[2], xend = ex_scenarios[2], 
                                 y = ex_projections[2], yend = true_obs),
               aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#579d1c") + 
  geom_text(data = data.frame(x = ex_scenarios[1]+0.03, 
                              y = mean(c(ex_projections[1], true_obs))), 
            aes(x = x, y = y), label = "observed\nbias", 
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

plot_grid(p_tot, p_scenario, p_miscal, nrow = 1)
ggsave("decompose_observed_error.pdf", width = 9, height = 3)

## multiple locations

obs_scenario_multloc = c(0.2, 0.45, 0.6, 0.75)
true_obs_multloc = max(true_projections) + diff(true_projections)*obs_scenario

ex_projections_multloc = list(
                      c(120, 70), 
                      c(130, 80), 
                      c(140, 60), 
                      c(110, 50)
                      )

true_projections_multloc = list(
  c(80, 60), 
  c(100, 70), 
  c(90, 40), 
  c(70, 30)
)

intermed_projection_multloc <- list()
true_obs_multloc <- list()

for(i in 1:length(obs_scenario_multloc)){
  intermed_projection_multloc[[i]] = max(ex_projections_multloc[[i]]) + diff(ex_projections_multloc[[i]])*obs_scenario_multloc[i]
  true_obs_multloc[[i]] = max(true_projections_multloc[[i]]) + diff(true_projections_multloc[[i]])*obs_scenario_multloc[i]
  example_data_tmp = data.frame(
    scenarios = rep(c(ex_scenarios, obs_scenario_multloc[i]),2), 
    projections = c(ex_projections_multloc[[i]], intermed_projection_multloc[[i]], true_projections_multloc[[i]], true_obs_multloc[[i]]),
    type = c("model", "model", "model", "observation", "observation", "observation"), 
    observable = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
  )
  ggplot(data = example_data_tmp) + 
    geom_vline(xintercept = ex_scenarios) + 
    geom_vline(xintercept = obs_scenario_multloc[i], linetype = "dotted") + 
    geom_line(aes(x = scenarios, y = projections, color = type), alpha = 0.3) + 
    # geom_segment(data = data.frame(x = obs_scenario, xend = obs_scenario, 
    #                                y = intermed_projection, yend = true_obs),
    #              aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#004586") + 
    # geom_segment(data = data.frame(x = ex_scenarios[1], xend = ex_scenarios[1],
    #                               y = ex_projections[1], yend = true_projections[1]),
    #             aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#83caff") +
    # geom_segment(data = data.frame(x = ex_scenarios[2], xend = ex_scenarios[2],
    #                                y = ex_projections[2], yend = true_projections[2]),
    #              aes(x = x, xend = xend, y = y, yend = yend), size = 2, color = "#83caff") +
    geom_point(aes(x = scenarios, y = projections, color = type, 
                   fill = interaction(type,observable)), shape = 21, size = 3) +
    guides(fill = FALSE) +
    scale_color_manual(values = c("black", "#ff420e"), 
                       labels = c("projection", "observation")) + 
    scale_fill_manual(values = c("white", "white", "black", "#ff420e")) + 
    scale_x_continuous(breaks = c(ex_scenarios, obs_scenario_multloc[i]), 
                       labels = c("low vax", "high vax", "realized"), 
                       limits = c(- 0.1, 1.1),
                       name = "scenario axis") + 
    scale_y_continuous(limits = range(list(true_projections_multloc, ex_projections_multloc)), 
                       name = "cumulative hospitalizations\n(projection axis)") + 
    theme_bw() + 
    theme(legend.position = "none", 
          panel.grid = element_blank())
  ggsave(paste0("scenario_error_components_l", i,".tiff"), width = 4, height = 3.5)
}

# now plot scenario axis vs. miscalibration 

miscal_multloc <- sapply(1:length(intermed_projection_multloc), function(i){intermed_projection_multloc[[i]] - true_obs_multloc[[i]]})


ribbon = predict(lm(miscal~scenario, data = data.frame(scenario = obs_scenario_multloc, 
                               miscal = miscal_multloc)), data.frame(scenario = ex_scenarios), se.fit = TRUE)
ribbon.df = data.frame(x = ex_scenarios, 
                       mean = ribbon$fit, 
                       lwr = ribbon$fit - 1.96*ribbon$se.fit, 
                       upr = ribbon$fit + 1.96*ribbon$se.fit)

ggplot() +
  geom_vline(xintercept = ex_scenarios) + 
  geom_hline(yintercept = 0) + 
  geom_ribbon(data = ribbon.df, 
            aes(x = x, ymin = lwr, ymax =upr), alpha = 0.2, fill = "#004586") + 
  geom_line(data = ribbon.df, 
            aes(x = x, y = mean), linetype = "dotted", color = "#004586") + 
  # geom_rect(data = data.frame(x = ex_scenarios[1],
  #                             xend = ex_scenarios[2], 
  #                               ymin = rep(min(miscal_multloc)*0.85, 2), 
  #                               ymax = rep(max(miscal_multloc)*1.15, 2)), 
  #             aes(xmin = x, xmax = xend, ymin = ymin, ymax = ymax), alpha = 0.2) + 
  geom_point(data = data.frame(scenario = obs_scenario_multloc, 
                               miscal = miscal_multloc), 
             aes(x = scenario, y = miscal), size = 3, color = "#004586") + 
  scale_x_continuous(breaks = c(ex_scenarios), 
                     labels = c("low vax", "high vax"), 
                     limits = c(-0.1, 1.1),
                     name = "scenario axis") + 
  scale_y_continuous(limits = c(-1,1)*max(ribbon.df),
                     name = "miscalibration error") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
ggsave("miscalregression.tiff", width = 3.5, height = 3)

