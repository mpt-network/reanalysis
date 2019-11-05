
library("MPTmultiverse")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
tmpe <- new.env()
source("fun_prep.R", local = tmpe)
check_set <- get(x = "check_set",envir = tmpe)
library("DescTools") # for CCC

load("all_pairs_core.RData")
source("fun_analysis.R")

##############

library("afex")
library("emmeans")
library("effects")
set_sum_contrasts()
library("sjPlot")
emm_options(lmer.df = "asymptotic")

##############

dv_trait <- all_pars_a4 %>% 
  filter(cond_x == "Trait PP")

dv_cmle <- all_pars_a4 %>% 
  filter(cond_x == "Comp MLE")

dv_cmle_nopc <- all_pars_a4 %>% 
  filter(cond_x == "Comp MLE") %>% 
  filter(!(model %in% "pc"))

###############

pd <- ggplot(dv_cmle, aes(x = x, y = adev)) +
  geom_point(alpha = 0.2) +
  facet_wrap("cond_y") +
  ylab("Absolute Deviation") +
  xlab("Complete Pooling MLE Estimate") + 
  theme_bw(base_size = 20)
ggsave(plot = pd, 
       filename = "figures/abs_dev_cmle_data.png", 
       width = 24, height = 19, units = "cm", 
       dpi = 500)

mm_cmle_polys <- fit_polynomials(
  formula = adev ~ poly(x, DEG)*cond_y + (1|condition:dataset), 
  data = dv_cmle, REML = FALSE, max_degree = 10)

do.call(anova, mm_cmle_polys)
# Data: dv_cmle
# Models:
# MODEL1: adev ~ poly(x, 1L) * cond_y + (1 | condition:dataset)
# MODEL2: adev ~ poly(x, 2L) * cond_y + (1 | condition:dataset)
# MODEL3: adev ~ poly(x, 3L) * cond_y + (1 | condition:dataset)
# MODEL4: adev ~ poly(x, 4L) * cond_y + (1 | condition:dataset)
# MODEL5: adev ~ poly(x, 5L) * cond_y + (1 | condition:dataset)
# MODEL6: adev ~ poly(x, 6L) * cond_y + (1 | condition:dataset)
# MODEL7: adev ~ poly(x, 7L) * cond_y + (1 | condition:dataset)
# MODEL8: adev ~ poly(x, 8L) * cond_y + (1 | condition:dataset)
# MODEL9: adev ~ poly(x, 9L) * cond_y + (1 | condition:dataset)
# MODEL10: adev ~ poly(x, 10L) * cond_y + (1 | condition:dataset)
#         Df    AIC    BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# MODEL1  18 -35308 -35173  17672   -35344                              
# MODEL2  26 -36998 -36804  18525   -37050 1706.36      8  < 2.2e-16 ***
# MODEL3  34 -37390 -37136  18729   -37458  407.76      8  < 2.2e-16 ***
# MODEL4  42 -38172 -37858  19128   -38256  797.58      8  < 2.2e-16 ***
# MODEL5  50 -38497 -38124  19299   -38597  341.88      8  < 2.2e-16 ***
# MODEL6  58 -39274 -38841  19695   -39390  792.24      8  < 2.2e-16 ***
# MODEL7  66 -39564 -39071  19848   -39696  306.02      8  < 2.2e-16 ***
# MODEL8  74 -39997 -39444  20072   -40145  449.07      8  < 2.2e-16 ***
# MODEL9  82 -40197 -39585  20181   -40361  216.38      8  < 2.2e-16 ***
# MODEL10 90 -40442 -39769  20311   -40622  260.32      8  < 2.2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

p_cmle_poly <- sjPlot::get_model_data(mm_cmle_polys[[10]], type = "pred", 
                                      terms = c("x[all]", "cond_y"))
p_cmle_poly$cond_y <- p_cmle_poly$group
pd + 
  geom_line(data = p_cmle_poly, mapping = aes(x = x, y = predicted), 
            color = "blue", size = 2)
ggsave(filename = "figures/abs_dev_cmle_data_poly.png", 
       width = 24, height = 19, units = "cm", 
       dpi = 500)

##---------------------------------------------------------------
##                    Only Latent-Trait as IV                   -
##---------------------------------------------------------------


dv_cmle_lt_nopc <- all_pars_a4 %>% 
  filter(cond_x == "Comp MLE") %>%
  filter(cond_y == "Trait PP") 

mm_cmle_lt_m <- lmer(adev ~ model + (1|orig_condition:dataset),
                     data = dv_cmle_lt_nopc, REML = FALSE)
mm_cmle_lt_p <- lmer(adev ~ parameter + (1|orig_condition:dataset),
                     data = dv_cmle_lt_nopc, REML = FALSE)
mm_cmle_lt_0 <- fit_polynomials(
  formula = adev ~ poly(x, DEG) + (1|orig_condition:dataset), 
  data = dv_cmle_lt_nopc, REML = FALSE, max_degree = 10)
mm_cmle_lt_1 <- fit_polynomials(
  formula = adev ~ poly(x, DEG) + parameter + (1|orig_condition:dataset), 
  data = dv_cmle_lt_nopc, REML = FALSE, max_degree = 10)
do.call(anova, c(mm_cmle_lt_m, mm_cmle_lt_p, mm_cmle_lt_0, mm_cmle_lt_1))
# Data: dv_cmle_lt_nopc
# Models:
# MODEL3: adev ~ poly(x, 1L) + (1 | orig_condition:dataset)
# MODEL4: adev ~ poly(x, 2L) + (1 | orig_condition:dataset)
# MODEL5: adev ~ poly(x, 3L) + (1 | orig_condition:dataset)
# MODEL6: adev ~ poly(x, 4L) + (1 | orig_condition:dataset)
# MODEL7: adev ~ poly(x, 5L) + (1 | orig_condition:dataset)
# MODEL8: adev ~ poly(x, 6L) + (1 | orig_condition:dataset)
# MODEL1: adev ~ model + (1 | orig_condition:dataset)
# MODEL9: adev ~ poly(x, 7L) + (1 | orig_condition:dataset)
# MODEL10: adev ~ poly(x, 8L) + (1 | orig_condition:dataset)
# MODEL11: adev ~ poly(x, 9L) + (1 | orig_condition:dataset)
# MODEL12: adev ~ poly(x, 10L) + (1 | orig_condition:dataset)
# MODEL2: adev ~ parameter + (1 | orig_condition:dataset)
# MODEL13: adev ~ poly(x, 1L) + parameter + (1 | orig_condition:dataset)
# MODEL14: adev ~ poly(x, 2L) + parameter + (1 | orig_condition:dataset)
# MODEL15: adev ~ poly(x, 3L) + parameter + (1 | orig_condition:dataset)
# MODEL16: adev ~ poly(x, 4L) + parameter + (1 | orig_condition:dataset)
# MODEL17: adev ~ poly(x, 5L) + parameter + (1 | orig_condition:dataset)
# MODEL18: adev ~ poly(x, 6L) + parameter + (1 | orig_condition:dataset)
# MODEL19: adev ~ poly(x, 7L) + parameter + (1 | orig_condition:dataset)
# MODEL20: adev ~ poly(x, 8L) + parameter + (1 | orig_condition:dataset)
# MODEL21: adev ~ poly(x, 9L) + parameter + (1 | orig_condition:dataset)
# MODEL22: adev ~ poly(x, 10L) + parameter + (1 | orig_condition:dataset)
#         Df     AIC     BIC logLik deviance    Chisq Chi Df Pr(>Chisq)    
# MODEL3   4 -4385.9 -4364.5 2197.0  -4393.9                               
# MODEL4   5 -4442.7 -4415.9 2226.3  -4452.7  58.7765      1  1.766e-14 ***
# MODEL5   6 -4520.8 -4488.7 2266.4  -4532.8  80.1319      1  < 2.2e-16 ***
# MODEL6   7 -4569.4 -4531.9 2291.7  -4583.4  50.5750      1  1.147e-12 ***
# MODEL7   8 -4579.2 -4536.4 2297.6  -4595.2  11.8271      1  0.0005838 ***
# MODEL8   9 -4663.3 -4615.1 2340.7  -4681.3  86.0705      1  < 2.2e-16 ***
# MODEL1  10 -4684.6 -4631.0 2352.3  -4704.6  23.2411      1  1.429e-06 ***
# MODEL9  10 -4704.3 -4650.7 2362.2  -4724.3  19.7695      0  < 2.2e-16 ***
# MODEL10 11 -4740.8 -4681.8 2381.4  -4762.8  38.4360      1  5.658e-10 ***
# MODEL11 12 -4767.4 -4703.1 2395.7  -4791.4  28.6208      1  8.803e-08 ***
# MODEL12 13 -4794.2 -4724.5 2410.1  -4820.2  28.7934      1  8.053e-08 ***
# MODEL2  50 -5170.2 -4902.3 2635.1  -5270.2 449.9962     37  < 2.2e-16 ***
# MODEL13 51 -5170.9 -4897.7 2636.4  -5272.9   2.7257      1  0.0987448 .  
# MODEL14 52 -5178.4 -4899.8 2641.2  -5282.4   9.5368      1  0.0020139 ** 
# MODEL15 53 -5224.8 -4940.8 2665.4  -5330.8  48.3514      1  3.563e-12 ***
# MODEL16 54 -5258.4 -4969.1 2683.2  -5366.4  35.5836      1  2.443e-09 ***
# MODEL17 55 -5274.3 -4979.6 2692.1  -5384.3  17.8927      1  2.337e-05 ***
# MODEL18 56 -5386.8 -5086.8 2749.4  -5498.8 114.5524      1  < 2.2e-16 ***
# MODEL19 57 -5430.9 -5125.5 2772.4  -5544.9  46.0953      1  1.126e-11 ***
# MODEL20 58 -5475.9 -5165.2 2796.0  -5591.9  47.0191      1  7.030e-12 ***
# MODEL21 59 -5508.3 -5192.2 2813.1  -5626.3  34.3502      1  4.604e-09 ***
# MODEL22 60 -5546.9 -5225.5 2833.5  -5666.9  40.6457      1  1.825e-10 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

map(c(mm_cmle_lt_p, mm_cmle_lt_0, mm_cmle_lt_1), performance::r2)

mm_cmle_lt_p_lm <- lm(adev ~ parameter, data = dv_cmle_lt_nopc)
mm_cmle_lt_0_lm <- fit_polynomials_lm(
  formula = adev ~ poly(x, DEG), 
  data = dv_cmle_lt_nopc, max_degree = 10)
mm_cmle_lt_1_lm <- fit_polynomials_lm(
  formula = adev ~ poly(x, DEG) + parameter, 
  data = dv_cmle_lt_nopc, max_degree = 10)

lm_sel <- tibble(
  model = map_chr(c(list(mm_cmle_lt_p_lm), mm_cmle_lt_0_lm, mm_cmle_lt_1_lm), 
              ~deparse(.$call$formula)), 
  aic = map_dbl(c(list(mm_cmle_lt_p_lm), mm_cmle_lt_0_lm, mm_cmle_lt_1_lm), AIC),
  bic = map_dbl(c(list(mm_cmle_lt_p_lm), mm_cmle_lt_0_lm, mm_cmle_lt_1_lm), BIC),
  r2 = map_dbl(c(list(mm_cmle_lt_p_lm), mm_cmle_lt_0_lm, mm_cmle_lt_1_lm), 
               ~summary(.)$r.squared)
)
lm_sel %>% 
  arrange(aic)

lm_sel %>% 
  arrange(bic)

lm_sel %>% 
  arrange(desc(r2)) %>% 
  print(n = Inf)


### IVs:

## "model"
ivs <- c("parameter", "log1p_hetero", "chisq", "")



##################################################################
##                       Regression Trees                       ##
##################################################################

### party package

library("partykit")
library("ggparty")

tree1 <- lmtree(adev ~ poly(x, 2) | 
                  log1p_hetero + chisq + rho + fungi + 
                  log1p_fit_x + log1p_fit_y +
                  p_fit_x + p_fit_y +
                  se_x + se_y + 
                  rel_par_weight_x + rel_par_weight_y +
                  rel_n_x + rel_n_y, data = dv_cmle_lt_nopc)
tree1

ggparty()
autoplot(tree1, show_fit = FALSE)


tree2 <- lmtree(adev ~ poly(x, 2) | 
                  log1p_hetero + chisq + rho + fungi + 
                  log1p_fit_x + log1p_fit_y +
                  p_fit_x + p_fit_y +
                  se_x + se_y + 
                  rel_par_weight_x + rel_par_weight_y +
                  rel_n_x + rel_n_y, data = dv_cmle_lt_nopc, 
                minsize = 100)
tree2

autoplot(tree2)

pred_df <- get_predictions(tree2, ids = "terminal", newdata =  function(x) {
  #browser()
  x
  # data.frame(
  #   citations = 1,
  #   price = exp(seq(from = min(x$`log(price/citations)`),
  #                   to = max(x$`log(price/citations)`),
  #                   length.out = 100)))
})

ggparty(tree2) +
  geom_edge() +
  geom_edge_label(max_length = 3) +
  geom_node_splitvar() +
  geom_node_plot(gglist = list(geom_point(aes(x = `poly(x, 2)`,
                                              y = adev),
                                          alpha = 0.8),
                               theme_bw(base_size = 10)))

tree3 <- lmtree(adev ~ poly(x, 2) | 
                  log1p_hetero + chisq + rho + fungi + 
                  log1p_fit_x + log1p_fit_y +
                  p_fit_x + p_fit_y +
                  se_x + se_y + 
                  rel_par_weight_x + rel_par_weight_y +
                  rel_n_x + rel_n_y, data = dv_cmle_lt_nopc, 
                minsize = 100)