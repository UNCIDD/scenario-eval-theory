library(ggplot2)
library(paletteer)
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
w = 0.25
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
b1 = 38
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


ggplot(data = plot_df, aes(x = scenario_value, y = projection_value, color = model)) + 
  geom_line(linewidth = 0.5, alpha = 0.7) + 
  geom_point(data = plot_df %>% filter(scenario_value %in% ex_scenarios, model == "truth"), size = 2, shape = 21, fill = "white") +
  geom_point(data = plot_df %>% filter(scenario_value %in% ex_scenarios, model != "truth"), size = 2) + 
  geom_point(data = plot_df %>% filter(scenario_value == obs_scenario, model == "truth"), size = 2) + 
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual(values = c(paletteer_d("nbapalettes::hornets2"), "red")) +
  # scale_colour_paletteer_d("wesanderson::Darjeeling1", direction = -1) +
  scale_x_continuous(breaks = c(ex_scenarios, obs_scenario), 
                     labels = c("low", "high", "realized"), 
                     name = "vaccine uptake (scenario axis)") + 
  scale_y_continuous(name = "final epidemic size\n(projection axis)") + 
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

ggsave("R/decision-theory/figures/projection_dv_examples.pdf", width = 3.5, height = 3.75)
  
