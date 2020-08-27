
library("afex")
library("DescTools")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))
source("fun_analysis_paper.R")
INCLUDE_GAM <- TRUE

load("all_pairs_core.RData")
all_pairs <- all_pairs %>% 
  mutate(rel_weight = (rel_par_weight_x + rel_par_weight_y)/2,
         rel_n = (rel_n_x + rel_n_y) / 2) %>% 
  mutate(se_x_w = if_else(se_x > 0.25, 0.25, se_x),
         se_y_w = if_else(se_y > 0.25, 0.25, se_y),
         rel_n_w = if_else(rel_n > 15000, 15000, rel_n)) %>% 
    mutate(cond_x2 = cond_y,
         cond_y2 = cond_x)


##----------------------------------------------------------------
##                          Pairs Plot                           -
##----------------------------------------------------------------

plot_text <- all_pairs %>%
  group_by(cond_x2, cond_y2) %>% 
  summarise(ccc = substr(formatC(
    CCC(x, y, na.rm = TRUE)$rho.c$est, 
    digits = 2, format = "f"), 2, 4))
plot_text$ccc

ppairs <- all_pairs %>% 
  #filter(!(model %in% "pm")) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.2) + #aes(size = trials)
  facet_grid(cond_x2~ cond_y2, switch = "both", as.table = TRUE) +
  geom_text(data=plot_text,
            aes(x = 0.15, y = 0.9, label=ccc),
            parse = FALSE, inherit.aes=FALSE, size = 5) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25), 
                     labels = c("0", "", "0.5", "", "1")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                     labels = c("0  ", "", "0.5", "", "1  ")) +
  #scale_size(range = c(0.5, 2.5)) +
  geom_smooth(method = "gam", se = FALSE, 
                formula = y ~ s(x, bs = "ts"), colour = "red") +
  labs(x = "", y = "") +
  theme(legend.position = "none") + 
  theme(panel.grid.major.x = element_line(), 
        panel.grid.minor = element_blank())
ggsave("figures_man/pairsplot_all.png", plot = ppairs,
       width = 29, height = 28, units = "cm", 
       dpi = 500)


##----------------------------------------------------------------
##                    Univariate Relationships                   -
##----------------------------------------------------------------

theme_set(theme_bw(base_size = 13) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))

sel_methods <- c("Comp MLE", "No asy", "Beta PP", "Trait PP")
#, "Trait_u PP"

## see corresponding Report in docs folder for full overview
str(all_pairs)
targ_cmle <- all_pairs %>% 
  filter(cond_y == "Comp MLE") %>% 
  filter(cond_x != "Comp MLE")
targ_lpp <- all_pairs %>% 
  filter(cond_y == "Trait PP") %>% 
  filter(cond_x != "Trait PP")

targ_both <- bind_rows(targ_cmle, targ_lpp) %>% 
  filter(cond_x %in% sel_methods) %>% 
  droplevels %>% 
  mutate(cond_label = factor(paste("DV:", cond_y), 
                             levels = c("DV: Comp MLE", "DV: Trait PP"))) %>% 
  mutate(cond_iv_label = factor(paste("IV:", cond_x), levels = 
                                  paste(paste("IV:", levels(cond_x)))))
  
ylab <- "Abs. deviation"


### parameter-level covariates
pest <- compare_continuous_covariate(data = targ_both, covariate = poly(y, 2), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Estimate (quadratic)")
psex <- compare_continuous_covariate(data = targ_both, covariate = se_x_w, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (IV)")
psey <- compare_continuous_covariate(data = targ_both, covariate = se_y_w, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (DV)")
psd <- compare_continuous_covariate(data = targ_both, covariate = sd_emp_inv, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Individual variability (SD)")
prho <- compare_continuous_covariate(data = targ_both, covariate = rho_med, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter correlations (median)")
pfungi <- compare_continuous_covariate(data = targ_both, covariate = fungi_max, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter trade-off (max)")
prelw <- compare_continuous_covariate(data = targ_both, covariate = log(rel_weight), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Relative information (log)")
preln <- compare_continuous_covariate(data = targ_both, covariate = log(rel_n_w), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Relative N (log)")

## data-set-level covariates
phetero <- compare_continuous_covariate(data = targ_both, covariate = log1p(p_hetero), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Hetereogeneity (log p + 1)")
pfit <- compare_continuous_covariate(data = targ_both, covariate = log1p(p_fit_x), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Model fit (IV, log p + 1)")


puniv_good <- cowplot::plot_grid(
  psey, 
  psex + theme(strip.text = element_blank()), 
  pfungi + theme(strip.text = element_blank()), 
  prelw + theme(strip.text = element_blank()), 
  preln + theme(strip.text = element_blank()), 
  ncol = 1, rel_heights = c(1.3, rep(1, 4)))
ggsave("figures_man/univariate_good.png", plot = puniv_good, 
       width = 28, height = 25, units = "cm", 
       dpi = 500)

puniv_notgood <- cowplot::plot_grid(
  pest, 
  psd + theme(strip.text = element_blank()), 
  prho + theme(strip.text = element_blank()), 
  phetero + theme(strip.text = element_blank()), 
  pfit + theme(strip.text = element_blank()), 
  ncol = 1, rel_heights = c(1.3, rep(1, 4)))
ggsave("figures_man/univariate_notgood.png", plot = puniv_notgood, 
       width = 28, height = 25, units = "cm", 
       dpi = 500)


##---------------------------------------------------------------
##                      Table of Covariates                     -
##---------------------------------------------------------------

sing_dat <- all_pairs %>% 
  filter(cond_y == "Comp MLE") %>% 
  filter(cond_x == "Trait PP")

ns <- sing_dat %>% 
  group_by(model2) %>% 
  summarise(n = n())

meansd <- sing_dat %>% 
  group_by(model2) %>% 
  summarise(across(c(sd_emp_inv, rho_med, fungi_max, 
                     rel_weight, rel_n_w, p_hetero, 
                     p_fit_x), .fns = ~ paste0(
                       formatC(mean(.), format = "f", digits = 2),
                       " ±",
                       formatC(sd(.), format = "f", digits = 2)
                     )))

minmax <- sing_dat %>% 
  group_by(model2) %>% 
  summarise(across(c(sd_emp_inv, rho_med, fungi_max, 
                     rel_weight, rel_n_w, p_hetero, 
                     p_fit_x), .fns = ~ paste0(
                       formatC(min(.), format = "f", digits = 2),
                       " - ", 
                       formatC(max(.), format = "f", digits = 2)
                     )))

left_join(ns, meansd)

left_join(ns, minmax)

##----------------------------------------------------------------
##              Correlation of Relevant Covariates               -
##----------------------------------------------------------------




rel_vars <- c("se_x_w", "se_y_w", "fungi_max", "rel_weight", "rel_n")

cors <- targ_cmle %>% 
  group_by(cond_x)  %>% 
  mutate(rel_weight = log(rel_weight),
         rel_n = log(rel_n)) %>% 
  select(rel_vars) %>%
  nest() %>% 
  summarise(cor = map(data, 
                      ~as_tibble(cor(., use = "na.or"), rownames = "par1"))) %>% 
  mutate(cor = map(cor, ~pivot_longer(., -par1, names_to = "par2", values_to = "cor"))) %>% 
  unnest(cor) %>% 
  filter(par1 != par2) %>% 
  pivot_wider(c(par1, par2), names_from = "cond_x", values_from = "cor") 

## filter out all repeating columns
cors %>% 
  filter(par2 != "se_x", par1 != "rel_n") %>% 
  slice(-8, -11, -12) 
  
