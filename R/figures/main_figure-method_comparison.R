library(ggplot2)
library(dplyr)

#### SETUP ---------------------------------------------------------------------
source("R/simulation/0-parameters.R")

all_ests = readRDS("output/simulation/estimated_errors_across_locations.rds")
method_comparsion_results = readRDS("output/simulation/performance_evaluation_results_across_locations.rds")

#### PANEL A: ALL ESTIMATED DISTRIBUTIONS --------------------------------------
approach_general_labs = c("approach 1: most plausible scenario", 
                          "approach 2: estimate observations", 
                          "approach 3: estimate error", "truth")
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
  scale_color_manual(values = c(RColorBrewer::brewer.pal(6, "Paired")[c(2, 4, 3, 6, 5)], "darkgray"), 
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

#### PANEL B: PERFORMANCE RESULTS ----------------------------------------------

method_comparsion_results = method_comparsion_results %>%
  left_join(all_ests_w_rank %>% select(model_id, scenario_id, approach, rank, true_rank, correct_flag))

p3 = ggplot(data = method_comparsion_results, aes(x = reorder(model_id, true_rank), y = mae, color = approach)) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_grid(rows = vars(scenario_id),
             labeller = labeller(scenario_id = scenario_labs), switch = "y") +
  labs(x = "model", y = "absolute difference between means") +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(6, "Paired")[c(2, 4, 3, 6, 5)], "darkgray")) +
  theme_bw(base_size = 7) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")
p3

p4 = ggplot(data = method_comparsion_results, aes(x = reorder(model_id, true_rank), y = ks_test_stat, color = approach)) + 
  geom_hline(data = method_comparsion_results %>%
               dplyr::select(n_df, n_truth) %>%
               unique() %>%
               mutate(ks_sig = 1.358*sqrt((n_df + n_truth)/(n_df*n_truth))), # 5% p value level 
             aes(yintercept = ks_sig), linewidth = 0.3) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  facet_grid(rows = vars(scenario_id),
             labeller = labeller(scenario_id = scenario_labs), switch = "y") + 
  labs(x = "model", y = "Kolmogorov-Smirnov test statistic") +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(6, "Paired")[c(2, 4, 3, 6, 5)], "darkgray")) +
  theme_bw(base_size = 7) + 
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside")
p4 

#### COMBINE INTO SINGLE PANEL -------------------------------------------------
# l = cowplot::get_legend(p1)
# cowplot::plot_grid(
#   p1 + theme(legend.position = "none"), 
#   cowplot::plot_grid(p3, p4, nrow = 1), 
#   l,
#   ncol = 1, rel_heights = c(0.46, 0.46, 0.08), labels = c("A", "B", NA), label_size = 10
# )

cowplot::plot_grid(
  p1, 
  cowplot::plot_grid(p3, p4, nrow = 1), 
  ncol = 1, rel_heights = c(0.6, 0.40), labels = c("A", "B"), label_size = 10
)

ggsave("output/figures/method_comparison.pdf", width = 7, height = 6)






