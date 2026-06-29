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

sigmoid_utility = function(y, l, r, y0, c){
  return(l/(1+exp(-r*(y+y0))) + c)
}

colors = c("#D55E00","#56B4E9", "#009E73")
names(colors) = c("linear", "quadratic", "sigmoid")

j_val = rep(3.5*1e2,2) #(population in thousands)
k_val = rep(5*1e2,2) #(population in thousands)
l_val = rep(2.25e2,2)
r_val = rep(10,2)
y0_val = rep(-0.3,2) # sigmoid inflection at y = 0.3 (within the projection range)
# S1 = more aggressive (and more expensive) action than S2; action cost differs by 500
c_val = c(2.5*1e2, 2*1e2)

utility_fns = expand.grid(y = seq(0, 0.7, length.out = 100), scenario_id = c("S1", "S2")) %>%
  mutate(scenario_int = as.integer(substr(scenario_id,2,2))) %>%
  mutate(linear = linear_utility(y, j_val[scenario_int], c_val[scenario_int]), 
         quadratic = quadratic_utility(y, k_val[scenario_int], c_val[scenario_int]), 
         sigmoid = sigmoid_utility(y, l_val[scenario_int], r_val[scenario_int], 
                                           y0_val[scenario_int], c_val[scenario_int])) %>%
  dplyr::select(-scenario_int) %>%
  melt(c("y", "scenario_id"), variable.name = "utility_fn", value.name = "utility")

# plot 1: utility functions
p1 = ggplot(utility_fns %>% filter(scenario_id == "S1"),
            aes(x = y, y = utility, color = utility_fn)) +
  geom_line(linewidth = 0.7) +
  labs(x = "projected outcome", color = "utility function") +
  scale_color_manual(values = colors) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic(base_size = 7) +
  theme(legend.position = "bottom", 
        axis.line = element_line(linewidth = 0.25))
p1

#### 2: PROJECTIONS ------------------------------------------------------------
## Projections chosen so that BOTH cross-model and within-model disagreement appears:
##   M1 (severe epidemic under both actions) — linear + quadratic prefer S1, sigmoid prefers S2
##   M2 (mild epidemic under both actions)   — all three utilities prefer S2
##   T  (truth: S1 effective, S2 moderate)   — linear + sigmoid prefer S1, quadratic prefers S2
## So: different MODELS disagree (each utility flips its decision across models), and within
##     M1 and T different UTILITY FUNCTIONS disagree, with a different utility being the outlier
##     in each case (sigmoid in M1 — saturated, doesn't reward S1's reduction; quadratic in T —
##     S2's outcome already small enough that the quadratic penalty doesn't justify S1's cost).
projected_values = expand_grid(model_id = c("M1", "M2", "T"),
                               scenario_id = c("S1", "S2")) %>%
  mutate(mu = c(0.40, 0.55,   # M1: S1, S2
                0.15, 0.22,   # M2: S1, S2
                0.14, 0.34),  # T:  S1, S2
         sigma = c(0.001, 0.001, 0.0015, 0.0015, 0, 0)) %>%
  reframe(final_size_samp = as.vector(mvrnorm(ifelse(model_id == "T", 1, samp_size), mu, sigma)), 
          draw_id = 1:ifelse(model_id == "T", 1, samp_size), .by = c("scenario_id", "model_id")) 

projected_values_wide = projected_values %>%
  dcast(model_id + draw_id ~ scenario_id, value.var = "final_size_samp")

# ggplot(projected_values_wide, aes(x = S1, y = S2, color = model_id)) + 
#   geom_point() + 
#   # scale_x_continuous(limits = c(0,1)) +
#   # scale_y_continuous(limits = c(0,1)) + 
#   theme_bw()


#### 3: UTILITIES --------------------------------------------------------------
utilities = projected_values %>% 
  mutate(scenario_int = as.integer(substr(scenario_id,2,2))) %>%
  mutate(linear = as.vector(linear_utility(final_size_samp, j_val[scenario_int], c_val[scenario_int])), 
         quadratic = as.vector(quadratic_utility(final_size_samp, k_val[scenario_int], c_val[scenario_int])), 
         sigmoid = as.vector(sigmoid_utility(final_size_samp, l_val[scenario_int], r_val[scenario_int], y0_val[scenario_int], c_val[scenario_int]))
  ) %>%
  dplyr::select(-scenario_int) %>%
  melt(c("scenario_id", "model_id", "draw_id", "final_size_samp"), variable.name = "utility_fn", value.name = "utility") 

utilities_wide = utilities %>%
  dcast(model_id + draw_id + utility_fn ~ scenario_id, value.var = "utility")

p2 = ggplot(utilities_wide %>% filter(model_id != "T"), aes(x = S1, y = S2)) + 
  # geom_point(aes(color = utility_fn, fill = utility_fn), alpha = 0.2, size = 0.5) +
  # geom_point(data = utilities_wide %>% filter(model_id == "T"), size = 4) +
  # geom_text(data = utilities_wide %>% filter(model_id == "T"), size = 4, 
  #           label = "*", color = "white", vjust = 0.75) +
  stat_ellipse(data = utilities_wide %>% filter(model_id != "T"),
                 aes(fill = utility_fn), geom = "polygon", alpha = 0.2) +
  geom_abline() +
  facet_wrap(vars(paste0("model ", substr(model_id,2,2)))) +
  labs(x = "projected utility (action 1)", y = "projected utility (action 2)", color = "utility function") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  coord_equal() +
  theme_bw(base_size = 7) + 
  theme(axis.line = element_line(linewidth = 0.25), 
        legend.position = "none", 
        panel.grid = element_blank(), 
        strip.background = element_blank())
p2

#### 4: DECISION LINES IN PROJECTION SPACE -------------------------------------
get_decision_line_sig = function(y2, y1, r_val, l_val, y0_val, c_val){
  sigmoid_utility(y2, l = l_val[2], r = r_val[2], y0 = y0_val[2], c = c_val[2]) - 
    sigmoid_utility(y1, l = l_val[1], r = r_val[1], y0 = y0_val[1], c = c_val[1])
}

## Decision line in (y1, y2) projection space: solve u(y2; S2) = u(y1; S1) for y2.
decision_lines = data.frame(y1 = seq(0.001, 0.55, length.out = 200))
decision_lines = decision_lines %>%
  mutate(linear    = (j_val[1] * y1   + c_val[1] - c_val[2]) / j_val[2],
         quadratic = sqrt((k_val[1] * y1^2 + c_val[1] - c_val[2]) / k_val[2])) %>%
  mutate(sigmoid = tryCatch(
    uniroot(get_decision_line_sig, interval = c(-1e6,1e6),
            y1 = y1, r_val, l_val, y0_val, c_val)$root,
    error = function(e) NA_real_), .by = c("y1")) %>%
  melt(c("y1"), variable.name = "utility_fn", value.name = "y2")

p3 = ggplot(projected_values_wide, aes(x = S1, y = S2)) +
  stat_ellipse(data = projected_values_wide %>% filter(model_id != "T"),
                 aes(group = model_id), geom = "polygon", alpha = 0.2) +
  # geom_point(alpha = 0.2, size = 0.5) +
  geom_line(data = decision_lines, aes(x = y1, y = y2, color = utility_fn), linewidth = 0.7) +
  geom_text(data = data.frame(model = paste("model", 1:2), 
                              x =  projected_values_wide %>% 
                                filter(model_id != "T") %>%
                                summarize(m = quantile(S1, 0.98), .by = "model_id") %>%
                                pull(m),
                              y = projected_values_wide %>% 
                                filter(model_id != "T") %>%
                                summarize(m = quantile(S2, 0.02), .by = "model_id") %>%
                                pull(m) 
                              ), 
            aes(x = x, y = y, label = model)) + 
  labs(x = "projected outcome (action 1)", y = "projected outcome (action 2)", color = "utility function") +
  scale_color_manual(values = colors) +
  scale_x_continuous(expand = c(0,0), limit = c(0,0.7)) +
  scale_y_continuous(expand = c(0,0), limit = c(0,0.7)) +
  theme_bw(base_size = 7) + 
  theme(axis.line = element_line(linewidth = 0.25), 
        legend.position = "none", 
        panel.grid = element_blank())

plot_grid(
  plot_grid(p1, p2, ncol = 1, rel_heights = c(0.52, 0.48)),
  p3, nrow = 1)

ggsave("R/decision-theory/figures/projection_space.pdf", width = 6, heigh = 4)  

