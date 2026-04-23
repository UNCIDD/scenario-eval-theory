library(ggplot2)
library(paletteer)
library(RColorBrewer)
library(cowplot)
library(dplyr)

### set up simple example
# first set scenario values and observations
ex_scenarios = c(0, 1)
obs_scenario = 0.7
full_scenario_vals = unique(sort(c(seq(ex_scenarios[1], ex_scenarios[2], length.out = 25), obs_scenario)))
vec_len = length(full_scenario_vals)

# now set true values
v = 140
w = 0.1
true_obs = v*w^obs_scenario - 3
true_projection = v*w^full_scenario_vals - 3
true_values = true_projection[c(1, vec_len)]

# projection 1 = exact shape, wrong values
projection_1 = v*w^full_scenario_vals - 25

# projection 2 = exact projection values, very wrong shape
a = true_values[2] - true_values[1]
c = true_values[2] - a + 2 # + 2 to add a little extra to be off slightly
projection_2 = a*full_scenario_vals^2 + c

# projection 3 = exact relationship between x and y, but wrong shape
m = true_values[2] - true_values[1]
b = true_values[1] + 15 # add a little to shift it
projection_3 = m*full_scenario_vals + b

# projection 4 = observation exactly, projected values very wrong
b1 = 11
m1 = (true_obs - b1)/obs_scenario
projection_4 = m1*full_scenario_vals + b1
  
plot_df = data.frame(
  model = c(rep("truth", vec_len), 
            rep("projection 1", vec_len), 
            rep("projection 2", vec_len), 
            rep("projection 3", vec_len), 
            rep("projection 4", vec_len)), 
  scenario_value = rep(full_scenario_vals, 5),
  projection_value = c(true_projection, projection_1, projection_2, projection_3, projection_4)
)

model_colors =c ("#E69F00", "#56B4E9", "#009E73", "#CC79A7") #"#f0e442"

p1 = ggplot(data = plot_df, aes(x = scenario_value, y = projection_value, 
                           color = model, linewidth = model, linetype = model)) + 
  geom_line(alpha = 0.7) + 
  geom_point(data = plot_df %>% filter(scenario_value %in% ex_scenarios, model != "truth"), size = 2.4) + 
  geom_point(data = plot_df %>% filter(scenario_value %in% ex_scenarios, model == "truth"), size = 2.4, shape = 21, fill = "white") +
  geom_point(data = plot_df %>% filter(scenario_value == obs_scenario, model == "truth"), size = 2.4) + 
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual(values = c(model_colors, "black")) +
  scale_linetype_manual(values = c(rep("solid", 4), "longdash")) +
  scale_linewidth_manual(values = c(rep(0.6,4), 1.1)) +
  # scale_colour_paletteer_d("wesanderson::Darjeeling1", direction = -1) +
  scale_x_continuous(breaks = c(ex_scenarios, obs_scenario), 
                     labels = c("low", "high", "realized"), 
                     name = "vaccine uptake (scenario axis)") + 
  scale_y_continuous(name = "final epidemic size\n(projection axis)") + 
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave("R/decision-theory/figures/projection_dv_examples.pdf", p1, width = 3.5, height = 3.75)
  

table_df = data.frame(
  projection_name = paste0("projection ", 1:4), 
  projection = 1:4,
  direction = c(TRUE, TRUE, TRUE, FALSE),
  shape = c(TRUE, FALSE, FALSE, FALSE), 
  scenario_value = c(FALSE, TRUE, FALSE, FALSE), 
  relative_value = c(TRUE, FALSE, TRUE, FALSE), 
  forecast = c(FALSE, FALSE, FALSE, TRUE)
) %>% 
  melt(c("projection")) %>%
  mutate(variable = factor(variable, levels = rev(c("projection_name", "direction", "scenario_value", "relative_value", "shape", "forecast"))))

p2 = ggplot(data = table_df, aes(x = projection, y = variable, fill = value))+  #, y = variable, 
  geom_tile(alpha = 0.2) + 
  geom_text(data = table_df %>% filter(variable == "projection_name"), 
            aes(label = value, color = value), size = 2) + 
  geom_vline(data = data.frame(x = 0:5 + 0.5), aes(xintercept = x)) +
  geom_hline(data = data.frame(y = 0:6 + 0.5), aes(yintercept = y)) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = c("white", "white", "white", "white", "white", "black")) +
  # scale_x_discrete(position = "top") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), labels = rev(c("", "direction", "absolute\nvalue", "relative\nvalue", "shape", "forecast\naccuracy"))) +
  theme_bw(base_size = 7) + 
  theme(axis.text.y = element_text(size = 6),
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        plot.margin = margin(0.3, 0.8, 0.3, 0.8, "cm"))
  

plot_grid(p1 + theme(legend.position = "none"), p2, ncol = 1, 
          rel_heights = c(0.6, 0.4))
ggsave("R/decision-theory/figures/projection_dv_examples_wtable.pdf", width = 3.25, height = 4)

