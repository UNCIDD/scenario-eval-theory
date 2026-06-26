library(dplyr)
library(reshape2)
library(ggplot2)
library(MASS) 
library(ggforce)
library(tidyr)
library(cowplot)

#### 1: UTILITY FUNCTIONS ------------------------------------------------------

samp_size = 500

linear_utility = function(y, j, c){
  return(y*j + c) 
}

quadratic_utility = function(y, k, c){
  return(k*y^2 + c)
}

logarithmic_utility = function(y, l, r, y0, c){
  return(l/(1+exp(-r*(y+y0))) + c)
}

j_val = rep(3.5*1e3,2) #(population in thousands)
k_val = rep(5*1e3,2) #(population in thousands)
l_val = rep(2.25e3,2)
r_val = rep(10,2)
y0_val = rep(-0.5,2)
# S1 = more aggressive (and more expensive) action than S2; action cost differs by 500
c_val = c(2.5*1e3, 2*1e3)

utility_fns = expand.grid(y = seq(0, 1, length.out = 100), scenario_id = c("S1", "S2")) %>%
  mutate(scenario_int = as.integer(substr(scenario_id,2,2))) %>%
  mutate(linear = linear_utility(y, j_val[scenario_int], c_val[scenario_int]), 
         quadratic = quadratic_utility(y, k_val[scenario_int], c_val[scenario_int]), 
         logarithmic = logarithmic_utility(y, l_val[scenario_int], r_val[scenario_int], 
                                           y0_val[scenario_int], c_val[scenario_int])) %>%
  dplyr::select(-scenario_int) %>%
  melt(c("y", "scenario_id"), variable.name = "utility_fn", value.name = "utility")

# plot 1: utility functions
p1 = ggplot(utility_fns, aes(x = y, y = utility, color = utility_fn, linetype = scenario_id)) + 
  geom_line() +
  theme_classic() + 
  theme()
p1

#### 2: PROJECTIONS ------------------------------------------------------------
## Projections chosen so that:
##   M1 (pessimistic about S2) — all three utilities prefer S1
##   M2 (optimistic about S2)  — all three utilities prefer S2
##   T  (truth, in between)    — linear prefers S1 but quadratic and sigmoid prefer S2
## i.e. different MODELS lead to different decisions (M1 vs M2), and
##      different UTILITY FUNCTIONS lead to different decisions (under T).
projected_values = expand_grid(model_id = c("M1", "M2", "T"),
                               scenario_id = c("S1", "S2")) %>%
  mutate(mu = c(0.13, 0.48,   # M1: S1, S2
                0.15, 0.22,   # M2: S1, S2
                0.14, 0.34),  # T:  S1, S2
         sigma = c(0.001, 0.001, 0.0015, 0.0015, 0.0005, 0.0005)) %>%
  reframe(final_size_samp = as.vector(mvrnorm(samp_size, mu, sigma)), 
          draw_id = 1:samp_size, .by = c("scenario_id", "model_id")) 

projected_values_wide = projected_values %>%
  dcast(model_id + draw_id ~ scenario_id, value.var = "final_size_samp")

ggplot(projected_values_wide, aes(x = S1, y = S2, color = model_id)) + 
  geom_point() + 
  # scale_x_continuous(limits = c(0,1)) +
  # scale_y_continuous(limits = c(0,1)) + 
  theme_bw()


#### 3: UTILITIES --------------------------------------------------------------
utilities = projected_values %>% 
  mutate(scenario_int = as.integer(substr(scenario_id,2,2))) %>%
  mutate(linear = as.vector(linear_utility(final_size_samp, j_val[scenario_int], c_val[scenario_int])), 
         quadratic = as.vector(quadratic_utility(final_size_samp, k_val[scenario_int], c_val[scenario_int])), 
         logarithmic = as.vector(logarithmic_utility(final_size_samp, l_val[scenario_int], r_val[scenario_int], y0_val[scenario_int], c_val[scenario_int]))
  ) %>%
  dplyr::select(-scenario_int) %>%
  melt(c("scenario_id", "model_id", "draw_id", "final_size_samp"), variable.name = "utility_fn", value.name = "utility") 

utilities_wide = utilities %>%
  dcast(model_id + draw_id + utility_fn ~ scenario_id, value.var = "utility")


p2 = ggplot(utilities_wide, aes(x = S1, y = S2, color = utility_fn, fill = utility_fn)) + 
  # geom_point() + 
  stat_ellipse(geom = "polygon", alpha = 0.3) +
  geom_abline() +
  facet_wrap(vars(model_id)) +
  theme_bw()

#### 4: DECISION LINES IN PROJECTION SPACE -------------------------------------
get_decision_line_log = function(y2, y1, r_val, l_val, y0_val, c_val){
  logarithmic_utility(y2, l = l_val[2], r = r_val[2], y0 = y0_val[2], c = c_val[2]) - 
    logarithmic_utility(y1, l = l_val[1], r = r_val[1], y0 = y0_val[1], c = c_val[1])
}

decision_lines = data.frame(y1 = seq(0.001, 0.5, length.out = 100))
decision_lines = decision_lines %>% 
  mutate(linear = 1/j_val[1] * (j_val[2] * y1 + c_val[2] - c_val[1]), 
         quadratic = sqrt(1/k_val[1] * (k_val[2] * y1^2 + c_val[2] - c_val[1]))) %>%
  mutate(logarithmic =uniroot(get_decision_line_log, interval = c(-1e6,1e6),
                              y1 = y1, r_val, l_val, y0_val, c_val)$root, .by = c ("y1")) %>%
  melt(c("y1"), variable.name = "utility_fn", value.name = "y2")

p3 = ggplot(projected_values_wide, aes(x = S1, y = S2)) + 
  # geom_point() + 
  stat_ellipse(geom = "polygon", alpha = 0.3) +
  geom_line(data = decision_lines, aes(x = y1, y = y2, color = utility_fn)) +
  facet_wrap(vars(model_id)) +
  # scale_x_continuous(limits = c(0,1)) +
  # scale_y_continuous(limits = c(0,1)) + 
  theme_bw()

plot_grid(p1, p2, p3, ncol = 1)

ggsave("R/decision-theory/figures/projection_space.pdf", width = 6, heigh = 4)  

