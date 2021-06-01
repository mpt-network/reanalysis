
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
         cond_y2 = cond_x) %>% 
  mutate(se_c = (se_x_w + se_y_w)/2) 
all_pairs <- all_pairs %>%
  mutate(cond_co = apply(cbind(as.character(cond_x), 
            as.character(cond_y)), 1, 
      function(x) paste(sort(x), collapse = " - ")))


##----------------------------------------------------------------
##                          Pairs Plot                           -
##----------------------------------------------------------------

plot_text <- all_pairs %>%
  group_by(cond_x2, cond_y2) %>% 
  summarise(ccc = substr(formatC(
    CCC(x, y, na.rm = TRUE)$rho.c$est, 
    digits = 2, format = "f"), 2, 4))
# plot_text$ccc

plot_rmse <- all_pairs %>%
  group_by(cond_x2, cond_y2) %>% 
  summarise(rmse = substr(formatC(
    sqrt(mean((x - y)^2)), 
    digits = 3, format = "f"), 2, 5))

ppairs <- all_pairs %>% 
  #filter(!(model %in% "pm")) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.2) + #aes(size = trials)
  facet_grid(cond_x2~ cond_y2, switch = "both", as.table = TRUE) +
  geom_text(data=plot_rmse,
            aes(y = 0.13, x = 0.85, label=rmse),
            parse = FALSE, inherit.aes=FALSE, size = 4) +
  # geom_text(data=plot_rmse,
  #           aes(y = 0.1, x = 0.85, label=rmse),
  #           parse = FALSE, inherit.aes=FALSE, size = 5, color = "blue") +
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


## ecdf
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
library("ggthemes")

all_pairs %>% 
  filter(cond_x != "No Bayes", cond_y != "No Bayes") %>% 
  ggplot(aes(abs_dev, color = cond_y2)) +
  stat_ecdf(geom = "step", 
            mapping = aes(y = after_stat(1 - y))) +
  facet_wrap("cond_x2", as.table = TRUE, nrow = 2) +
  coord_cartesian(xlim = c(0.00, 0.35), expand = FALSE) + 
  labs(y = expression(paste(
    italic("Pr"), "(abs. deviation > ", italic("x"), ")")),  
       x = "Abs. deviation") +
  guides(color = guide_legend(title = NULL)) + scale_colour_colorblind()
ggsave("figures_man/eccdf.png", 
       width = 15, height = 12.5, units = "cm", 
       dpi = 500)

theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))


##------------------------
##  W/o no pooling Bayes  
##------------------------

unique_pairs <- all_pairs %>% 
  distinct(model, model2, dataset, parameter, cond_co, 
           parameter_o, condition, orig_condition, .keep_all = TRUE) %>% 
  filter(cond_x != "No Bayes", cond_y != "No Bayes", cond_x != cond_y) 

unique_pairs %>% 
  summarise(
    l05 = mean(abs_dev > .05),
    l10 = mean(abs_dev > .10),
    l25 = mean(abs_dev > .25),
    l33 = mean(abs_dev > .33),
    l50 = mean(abs_dev > .50),
    l75 = mean(abs_dev > .75),
    l90 = mean(abs_dev > .90),
  )

unique_pairs %>% 
  group_by(cond_co) %>% 
  summarise(
    max = max(abs_dev), 
    which = model2[which.max(abs_dev)],
    par = parameter[which.max(abs_dev)],
    data = dataset[which.max(abs_dev)]
  ) %>% 
  arrange(max) %>% 
  print(n = Inf)

unique_pairs %>% 
  group_by(cond_co) %>% 
  summarise(max = max(abs_dev), 
            which = model2[which.max(abs_dev)]) %>% 
  count(which) %>% 
  mutate(prob = n / sum(n))

unique_pairs %>% 
  summarise(mean = mean(abs_dev > 0.75),
         sum = sum(abs_dev > 0.75))

unique_pairs %>% 
  filter(abs_dev > 0.75) %>% 
  count(dataset, model, cond_co)

##-------------------
##  Partial Pooling  
##-------------------

pp_methods <- c("Beta PP", "Trait_u PP", "Trait PP")

unique_pairs %>% 
  filter(cond_x %in% pp_methods, cond_y %in% pp_methods) %>% 
  arrange(desc(abs_dev)) %>% 
  select(model:abs_dev, cond_x, cond_y) 

unique_pairs %>% 
  filter(cond_x %in% pp_methods, cond_y %in% pp_methods) %>% 
  arrange(desc(abs_dev)) %>% 
  select(model:abs_dev, cond_x, cond_y) %>% 
  filter(abs_dev > .1) %>% 
  count(model)


##---------------------
##  Bootstrap Methods  
##---------------------

boot_methods <- c("No PB", "No NPB")

all_pairs %>% 
  filter(cond_x %in% boot_methods, cond_y %in% boot_methods) %>% 
  arrange(desc(abs_dev)) %>% 
  slice(seq(1, n(), by = 2))


##---------------
##  All Methods  
##---------------

comp_methods <- c("Comp MLE", "Comp Bayes")

largest <- all_pairs %>% 
  arrange(desc(abs_dev)) %>% 
  select(model, model2, parameter, abs_dev, cond_x, cond_y) %>% 
  filter(abs_dev > 0.75)

## are there any large divergences not including complete pooling
largest %>% 
  filter(!(cond_x %in% comp_methods) & !(cond_y %in% comp_methods))

## 
largest %>% 
  mutate(comp_method = 
           case_when(
             (cond_x %in% comp_methods) & (cond_y %in% comp_methods) ~ "Comp MLE",
             cond_x %in% comp_methods ~ as.character(cond_x), 
             TRUE ~ as.character(cond_y)),
         other_method = case_when(
             (cond_x %in% comp_methods) & (cond_y %in% comp_methods) ~ "Comp Bayes",
             cond_x %in% comp_methods ~ as.character(cond_y), 
             TRUE ~ as.character(cond_x))) %>% 
  select(-cond_x, -cond_y) %>% 
  unique %>% 
  print(n = Inf)

all_pairs %>% 
  group_by(cond_x, cond_y) %>% 
  summarise(max = max(abs_dev),
            model = model[which.max(abs_dev)]) %>% 
  ungroup %>% 
  arrange(max) %>% 
  print(n = Inf)



##################################################################
##        Selected Pairs of Target and Predicting Method        ##
##################################################################


theme_set(theme_bw(base_size = 13) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))

sel_methods <- c("Comp MLE", "No asy", "Beta PP", "Trait PP")

rel_vars <- c("se_x_w", "se_y_w", "fungi_max", "rel_weight", "rel_n_w")
other_vars <- c("y", "sd_emp_inv", "rho_med", "p_hetero", "p_fit_x")
all_vars <- c(rel_vars, other_vars)
form_main <- paste0("abs_dev ~ ", paste(all_vars, collapse = " + "))

all_pairs_red <- all_pairs %>% 
  mutate(filtervar = create_filter_from_formula(form_main)) %>% 
  filter(filtervar)

1 - nrow(all_pairs_red)/nrow(all_pairs)

## see corresponding Report in docs folder for full overview
str(all_pairs)
targ_cmle <- all_pairs_red %>% 
  filter(cond_y == "Comp MLE") %>% 
  filter(cond_x != "Comp MLE")
targ_lpp <- all_pairs_red %>% 
  filter(cond_y == "Trait PP") %>% 
  filter(cond_x != "Trait PP")

targ_both <- bind_rows(targ_cmle, targ_lpp) %>% 
  filter(cond_x %in% sel_methods) %>% 
  droplevels %>% 
  mutate(cond_label = factor(paste("R:", cond_y), 
                             levels = c("R: Comp MLE", "R: Trait PP"))) %>% 
  mutate(cond_iv_label = factor(paste("C:", cond_x), levels = 
                                  paste(paste("C:", levels(cond_x)))))

ylab <- "Abs. deviation"
  

##----------------------------------------------------------------
##                  Absolute Mean Deviation Plot                 -
##----------------------------------------------------------------

targ_both %>% 
  ggplot(aes(x = cond_x, y = abs_dev)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_violin(fill = "transparent", width = 1.1) +
  stat_summary(fun = mean, fun.max = mean, fun.min = mean, fatten = 0.9) +
  facet_wrap("cond_label", scales = "free_x") +
  labs(x = "Comparison method", y = ylab)
ggsave("figures_man/mad.png", 
       width = 16, height = 7, units = "cm", 
       dpi = 500)



##----------------------------------------------------------------
##                    Univariate Relationships                   -
##----------------------------------------------------------------





### parameter-level covariates
pest <- compare_continuous_covariate(data = targ_both, covariate = poly(y, 2), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Value of parameter estimate (reference method, quadratic)")
psex <- compare_continuous_covariate(data = targ_both, covariate = se_x_w, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (comparison method)")
psey <- compare_continuous_covariate(data = targ_both, covariate = se_y_w, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (reference method)")
psd <- compare_continuous_covariate(data = targ_both, covariate = sd_emp_inv, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Individual variability (SD)")
prho <- compare_continuous_covariate(data = targ_both, covariate = rho_med, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter correlations (median)")
pfungi <- compare_continuous_covariate(data = targ_both, covariate = fungi_max, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter trade-offs (max)")
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
  xlab("Model fit (comparison method, log p + 1)")


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


##----------------------------------------------------------------
##                          RMSE table                           -
##----------------------------------------------------------------

rmse_tab_1 <- targ_both %>% 
  group_by(cond_label, cond_iv_label) %>% 
  summarise(mean = mean(abs_dev), 
            sd = sd(abs_dev), 
            rmse = sqrt(mean((mean(abs_dev) - abs_dev)^2)), 
            sigma = sigma(lm(abs_dev ~ 1))) %>% 
  mutate(across(where(is.double), ~formatC(.x, digits = 3, format = "f")))
rmse_tab_1

rmse_tab_2 <- targ_both %>% 
  group_by(cond_label, cond_iv_label) %>% 
  summarise(
    model = sigma(lm(abs_dev ~ model)),
    model2 = sigma(lm(abs_dev ~ model2)), 
    parameter = sigma(lm(abs_dev ~ parameter)),
    dataset = sigma(lm(abs_dev ~ dataset)), 
    population = sigma(lm(abs_dev ~ population)), 
    sci_goal = sigma(lm(abs_dev ~ sci_goal)) 
  ) %>% 
  mutate(across(where(is.double), ~formatC(.x, digits = 3, format = "f")))
rmse_tab_2

rmse_tab_3 <- targ_both %>% 
  group_by(cond_label, cond_iv_label) %>% 
  summarise(
    model = length(unique(model)),
    model2 = length(unique(model2)), 
    parameter = length(unique(parameter)),
    dataset = length(unique(dataset)), 
    population = length(unique(population)),
    sci_goal = length(unique(sci_goal))
  ) %>% 
  mutate(across(where(is.double), ~formatC(.x, digits = 3, format = "f")))
rmse_tab_3

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
  

#################################################################
##                     Plots in Discussion                     ##
#################################################################

all_pairs_red <- all_pairs_red %>% 
  mutate(par = str_remove(parameter, ".+:")) %>% 
  mutate(par = str_remove(par, "_")) %>% 
  mutate(par2 = case_when(
    model == "quad" ~ paste0("scriptstyle(",
                            str_extract(par, "[[:upper:]]+"), 
                            "[", 
                            str_extract(par, "[[:lower:]]+"), 
                            "])"),
    model == "hb" & nchar(par) > 1 ~ paste0("italic(", 
                            substr(par, 1, nchar(par)-1), 
                            ")[", 
                            toupper(substr(par, nchar(par), nchar(par))), 
                            "]"),
    nchar(par) > 1 ~ paste0("italic(", 
                            substr(par, 1, nchar(par)-1), 
                            "[", 
                            substr(par, nchar(par), nchar(par)), 
                            "])"),
    TRUE ~ paste0("italic(", par, ")")
  ))
# %>% 
#   select(par, par2) %>% 
#   {table(.$par2)}

all_pairs_red$mlab <- factor(
  all_pairs_red$model2, 
  levels = c("2htsm_4", "2htsm_5d","2htsm_6e","c2ht6","c2ht8",
             "pc","pd_s","pd_e","pm","hb","rm","real","quad"),
  labels = c("2HTSM 4", "2HTSM 5d","2HTSM 6e","c2HT 6","c2HT 8",
             "PC","PD","PD_e","PM","HB","RM","ReAL","QUAD")
)

##----------------------------------------------------------------
##            Absolute Mean Deviation Plot (Model)               -
##----------------------------------------------------------------

targ_cmle_lpp<- all_pairs_red %>% 
  filter(cond_y == "Comp MLE") %>% 
  filter(cond_x == "Trait PP") %>%
  filter(model2 != "pd_e")

targ_cmle_lpp %>% 
  ggplot(aes(x = par2, y = abs_dev)) +
  geom_boxplot(width = 0.08, outlier.shape = NA) +
  geom_violin(fill = "transparent", width = 0.8) +
  stat_summary(fun = mean, fun.max = mean, fun.min = mean, fatten = 0.9) +
  facet_wrap("mlab", scales = "free_x", ncol=4) +
  labs(x = "Parameter", y = ylab)+
  theme(axis.text.x = element_text(angle = 0,size=10)) +
  scale_x_discrete(labels = ggplot2:::parse_safe) 
# guide = guide_axis(check.overlap = TRUE)

ggsave("figures_man/mad_model.png", 
       width = 16, height = 16, units = "cm", 
       dpi = 500)

##----------------------------------------------------------------
##                          Trade off Plot                       -
##----------------------------------------------------------------

targ_cmle_lpp %>% 
  ggplot(aes(x = par2, y = fungi_max)) +
  geom_boxplot(width = 0.08, outlier.shape = NA) +
  geom_violin(fill = "transparent", width = 0.8) +
  stat_summary(fun = mean, fun.max = mean, fun.min = mean, fatten = 0.9) +
  facet_wrap("mlab", scales = "free_x", ncol=4) +
  labs(x = "Parameter", y = "Max. trade-off")+
  theme(axis.text.x = element_text(angle = 0,size=10))  +
  scale_x_discrete(labels = ggplot2:::parse_safe)


ggsave("figures_man/trade_model.png", 
       width = 16, height = 16, units = "cm", 
       dpi = 500)

##----------------------------------------------------------------
##                           SE Plot                             -
##----------------------------------------------------------------

targ_cmle_lpp %>% 
  ggplot(aes(x = par2, y = se_c)) +
  geom_boxplot(width = 0.08, outlier.shape = NA) +
  geom_violin(fill = "transparent", width = 0.8) +
  stat_summary(fun = mean, fun.max = mean, fun.min = mean, fatten = 0.9) +
  facet_wrap("mlab", scales = "free_x", ncol=4) +
  labs(x = "Parameter", y = "Standard error") +
  theme(axis.text.x = element_text(angle = 0,size=10))  +
  scale_x_discrete(labels = ggplot2:::parse_safe)

ggsave("figures_man/se_model.png", 
       width = 16, height = 16, units = "cm", 
       dpi = 500)

