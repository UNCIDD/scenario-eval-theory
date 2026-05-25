library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)

source("R/simulation/0-parameters.R")
source("R/simulation/0-helper-functions.R")


### plot example of bad fit (approach 2) ----
rep_id_choice = 122

all_ests = readRDS("output/simulation/estimated_errors_across_locations_sens.rds") %>% 
  filter(rep_id == rep_id_choice)

approach_labs_wrap = c("truth", "approach 1", 
                       "approach 2\n(no covariates)", "approach 2\n(covariates)", 
                       "approach 3\n(no covariates)", "approach 3\n(covariates)")
names(approach_labs_wrap) = names(approach_labs)


approach_general_labs = c("approach 1: most plausible scenario", 
                          "approach 2: estimate error", 
                          "approach 3: estimate observations", "truth")
names(approach_general_labs) = c(as.character(1:3), "t")

all_ests_w_rank = all_ests %>%
  filter(round(quantile,4) %in% round(c(0.05,0.25,0.5,0.75,0.95),4)) %>%
  mutate(quantile = paste0("Q", quantile*100)) %>%
  dcast(model_id + scenario_id + vax_cov + approach ~ quantile, value.var = "value") %>%
  mutate(rank = rank(Q50), .by = c("approach", "scenario_id")) %>%
  left_join(all_ests %>%
              filter(round(quantile,4) %in% round(c(0.05,0.25,0.5,0.75,0.95),4), scenario_id == "S1") %>%
              mutate(quantile = paste0("Q", quantile*100)) %>%
              dcast(model_id + scenario_id + vax_cov + approach ~ quantile, value.var = "value") %>%
              mutate(rank = rank(Q50), .by = c("approach", "scenario_id")) %>% filter(approach == "truth") %>%
              rename(true_rank = rank) %>% dplyr::select(model_id, true_rank)  
  ) %>%
  mutate(correct_flag = ifelse(rank == true_rank, TRUE, FALSE)) %>%
  mutate(correct_flag = ifelse(approach == "truth", FALSE, correct_flag)) %>%
  mutate(approach_general = substr(approach, 1,1), covariates_used = ifelse(grepl("-cov", approach), TRUE, FALSE))


sep_amount = 0.1
p1 = all_ests_w_rank %>%
  mutate(model_id_ranked = reorder(model_id, true_rank)) %>%
  filter(approach != "truth") %>%
  mutate(x_vals = as.integer(as.factor(model_id_ranked)) + ifelse(approach == 1, sep_amount, ifelse(covariates_used, 0, 2*sep_amount))) %>%
  ggplot(aes(color = approach)) +
  geom_segment(aes(x = x_vals, xend = x_vals, y = Q5, yend = Q95), linewidth = 0.3, position = position_nudge(x = sep_amount)) + 
  geom_segment(aes(x = x_vals, xend = x_vals, y = Q25, yend = Q75), linewidth = 0.7, position = position_nudge(x = sep_amount)) +
  geom_point(aes(x = x_vals, y = Q50), size = 1, position = position_nudge(x = sep_amount)) +
  geom_segment(data = all_ests_w_rank %>%
                 mutate(model_id_ranked = reorder(model_id, true_rank)) %>%
                 filter(approach == "truth") %>% dplyr::select(-approach_general)%>%
                 mutate(x_vals = as.integer(as.factor(model_id_ranked)) + -sep_amount/2),
               aes(x = x_vals, xend = x_vals,
                   y = Q5, yend = Q95, color = "truth"), linewidth = 0.3, position = position_nudge(x = -sep_amount)) +
  geom_segment(data = all_ests_w_rank %>%
                 mutate(model_id_ranked = reorder(model_id, true_rank)) %>%
                 filter(approach == "truth") %>% dplyr::select(-approach_general) %>%
                 mutate(x_vals = as.integer(as.factor(model_id_ranked)) + -sep_amount/2),
               aes(x = x_vals, xend = x_vals,
                   y = Q25, yend = Q75, color = "truth"), linewidth = 0.7, position = position_nudge(x = -sep_amount)) +
  geom_point(data = all_ests_w_rank %>%
               mutate(model_id_ranked = reorder(model_id, true_rank)) %>%
               filter(approach == "truth") %>% dplyr::select(-approach_general) %>%
               mutate(x_vals = as.integer(as.factor(model_id_ranked)) + -sep_amount/2),
             aes(x = x_vals, y = Q50, color = "truth"), size = 1, position = position_nudge(x = -sep_amount)) +
  facet_grid(cols = vars(approach_general), rows = vars(scenario_id), 
             labeller = labeller(scenario_id = scenario_labs, approach_general = approach_general_labs), switch = "y") + 
  labs(x = "model", y = "distribution of errors across locations", color = "error estimation approach") +
  scale_color_manual(breaks=c("1", "truth", "2-covariates", "2-nocovariates", "3-covariates", "3-nocovariates"), 
                     values = c(RColorBrewer::brewer.pal(6, "Paired")[2],
                                "darkgray",
                                RColorBrewer::brewer.pal(6, "Paired")[c(4, 3, 6, 5)]), 
                     labels = approach_labs) +
  scale_x_continuous(breaks = 1:n_models, 
                     labels = levels(all_ests_w_rank %>% mutate(model_id_ranked = reorder(model_id, true_rank)) %>% pull(model_id_ranked))) + 
  scale_shape_manual(values = c(NA, 8)) +
  scale_y_continuous(limits = 0.25*c(-1,1)) +
  theme_bw(base_size = 7) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(), 
        strip.placement = "outside") 
p1

# now plot fits
approach2_nocov_errors_across_locs = readRDS( "output/simulation/estimated_errors_approach2_full_across_locs_sens.rds") %>% 
  filter(rep_id == rep_id_choice, quantile %in% c(0.05, 0.5, 0.95)) %>%
  mutate(quantile = paste0("Q", quantile*100)) %>%
  dcast(vax_cov + rep_id + model_id + scenario_id ~ quantile) %>% 
  mutate(model_id = factor(model_id, levels = paste0("M", 1:10)))

sim_out = readRDS("output/simulation/sim_out_sens.rds")
error_df = sim_out[[rep_id_choice]]$errors %>% 
  mutate(model_id = factor(model_id, levels = paste0("M", 1:10)))

p2 = ggplot(data = error_df %>% filter(scenario_id == "T"), aes(x = vax_cov)) + 
  geom_point(aes(y = error), size = 1) + 
  geom_ribbon(data = approach2_nocov_errors_across_locs, 
              aes(ymin = Q5, ymax = Q95), alpha = 0.2) + 
  geom_line(data = approach2_nocov_errors_across_locs, aes(y = Q50)) + 
  facet_wrap(vars(model_id), ncol = 5) + 
  scale_x_continuous(name =  "realized vaccine uptake", breaks = c(0.3, 0.4, 0.5)) +
  theme_bw(base_size = 7) + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank())
plot_grid(p1, p2, ncol = 1, rel_heights = c(0.55, 0.45), labels = c("A", "B"))

ggsave("output/figures/supp_poor-fit-example.pdf", width = 8, height = 7)
