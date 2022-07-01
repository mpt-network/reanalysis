library("afex")
library("DescTools")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))
source("fun_analysis_paper.R")
source("docs/fun_eai.R")
INCLUDE_GAM <- TRUE

library("mgcv")
library("gratia")
library("cowplot")

load("all_pairs_core.RData")
levels(all_pairs$cond_x) <- new_method_labels(levels(all_pairs$cond_x))
levels(all_pairs$cond_y) <- new_method_labels(levels(all_pairs$cond_y))

all_pairs <- all_pairs %>% 
  mutate(rel_weight = (rel_par_weight_x + rel_par_weight_y)/2,
         rel_n = (rel_n_x + rel_n_y) / 2) %>% 
  mutate(se_x_w = if_else(se_x > 0.25, 0.25, se_x),
         se_y_w = if_else(se_y > 0.25, 0.25, se_y),
         rel_n_w = if_else(rel_n > 15000, 15000, rel_n)) %>% 
    mutate(cond_x2 = cond_y,
         cond_y2 = cond_x) %>% 
  mutate(prop_ns = rowMeans(cbind(prop_ns_nopb, prop_ns_nonpb, 
                                  prop_ns_noasy, prop_ns_trait), 
                            na.rm = TRUE)) %>% 
  #mutate(se_c = (se_x_w + se_y_w)/2) %>%
  mutate(se_c = rowMeans(cbind(se_x_w, se_y_w), na.rm = TRUE)) %>% 
  mutate(SAI = factor(if_else(model == "pc", "Not SAI", "SAI"), 
                      levels = c("SAI", "Not SAI"))) %>% 
  mutate(z_abs_dev = abs_dev/se_c) %>% 
  filter(parameter != "rm:g")
all_pairs <- all_pairs %>%
  mutate(cond_co = apply(cbind(as.character(cond_x), 
            as.character(cond_y)), 1, 
      function(x) paste(sort(x), collapse = " - ")))

all_pairs_nopc <- all_pairs %>% 
  filter(model != "pc") %>% 
  droplevels %>% 
  mutate(cond_label = factor(paste("R:", cond_y))) %>% 
  mutate(cond_iv_label = factor(paste("C:", cond_x)))

sel_methods <- c("CP-MLE", "NP-MLE", "PP-B", "PP-LT-C")


##################################################################
##                    Magnitude of Deviation                    ##
##################################################################

##----------------------------------------------------------------
##                          Pairs Plot                           -
##----------------------------------------------------------------

plot_rmse <- all_pairs %>%
  group_by(cond_x2, cond_y2) %>% 
  summarise(rmse = substr(formatC(
    sqrt(mean((x - y)^2)), 
    digits = 3, format = "f"), 2, 5))

plot_rmse_nopc <- all_pairs_nopc %>%
  group_by(cond_x2, cond_y2) %>% 
  summarise(rmse = substr(formatC(
    sqrt(mean((x - y)^2)), 
    digits = 3, format = "f"), 2, 5))

ppairs <- all_pairs %>% 
  #filter(!(model %in% "pm")) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(mapping = aes(shape = SAI), alpha = 0.2) + #aes(size = trials)
  facet_grid(cond_x2~ cond_y2, switch = "both", as.table = TRUE) +
  geom_text(data=plot_rmse,
            aes(y = 0.13, x = 0.85, label=rmse),
            parse = FALSE, inherit.aes=FALSE, size = 4) +
  geom_text(data=plot_rmse_nopc,
            aes(y = 0.9, x = 0.16, label=rmse),
            parse = FALSE, inherit.aes=FALSE, size = 4) +
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
ggsave("figures_man/pairsplot_all_v2.png", plot = ppairs,
       width = 29, height = 28, units = "cm", 
       dpi = 500)

all_rmse <- all_pairs %>%
   filter(as.numeric(cond_x) > as.numeric(cond_y)) %>% ## only select upper diag
  group_by(cond_x2, cond_y2) %>% 
  summarise(rmse = sqrt(mean((x - y)^2)),
            rmse_nopc = sqrt(mean( ((x - y)^2)[model != "pc"] ))) %>% 
  ungroup %>% 
  mutate(diff = rmse - rmse_nopc)

### SAI versus non-SAI

all_rmse %>% 
  filter(rmse > rmse_nopc) %>% 
  summarise(n = n(),
            mean = mean(diff))

all_rmse %>% 
  filter(rmse < rmse_nopc)

all_rmse %>% 
  filter(rmse < rmse_nopc) %>% 
  summarise(n = n(),
            mean = mean(diff))

#### RMSE overview
all_rmse %>% 
  summarise(across(
    .cols = c(rmse, rmse_nopc), 
    .fns = c("min" = min, "max" = max)
    ))

##---------------------------------------------------------------
##                          ECDF Plot                           -
##---------------------------------------------------------------

theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
library("ggthemes")

all_pairs_nopc %>% 
  filter(cond_x != "NP-Bayes", cond_y != "NP-Bayes") %>% 
  ggplot(aes(abs_dev, color = cond_y2)) +
  stat_ecdf(geom = "step", 
            mapping = aes(y = after_stat(1 - y))) +
  facet_wrap("cond_x2", as.table = TRUE, nrow = 2) +
  coord_cartesian(xlim = c(0.00, 0.35), expand = FALSE) + 
  labs(y = expression(paste(
    italic("Pr"), "(abs. deviation > ", italic("x"), ")")),  
       x = "Abs. deviation") +
  guides(color = guide_legend(title = NULL)) + scale_colour_colorblind() +
  scale_y_continuous()
ggsave("figures_man/eccdf.png", 
       width = 15, height = 12.5, units = "cm", 
       dpi = 500)

theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))

### eccdf analysis

unique_pairs_nopc <- all_pairs_nopc %>% 
  distinct(model, model2, dataset, parameter, cond_co, 
           parameter_o, condition, orig_condition, .keep_all = TRUE) %>% 
  filter(cond_x != "NP-Bayes", cond_y != "NP-Bayes", cond_x != cond_y) 

unique_pairs_nopc %>% 
  summarise(
    l05 = mean(abs_dev > .05),
    l10 = mean(abs_dev > .10),
    l25 = mean(abs_dev > .25),
    l33 = mean(abs_dev > .33),
    l50 = mean(abs_dev > .50),
    l75 = mean(abs_dev > .75),
    l90 = mean(abs_dev > .90),
  )

unique_pairs <- all_pairs %>% 
  distinct(model, model2, dataset, parameter, cond_co, 
           parameter_o, condition, orig_condition, .keep_all = TRUE) %>% 
   filter(cond_x != "NP-Bayes", cond_y != "NP-Bayes", cond_x != cond_y)

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

##---------------
##  Worst Cases  
##---------------

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

## worst cases without PC

unique_pairs_nopc %>% 
  group_by(cond_co) %>% 
  summarise(
    max = max(abs_dev), 
    which = model2[which.max(abs_dev)],
    par = parameter[which.max(abs_dev)],
    data = dataset[which.max(abs_dev)]
  ) %>% 
  arrange(max) %>% 
  print(n = Inf)


##################################################################
##           Testing Theoretically Derived Prediction           ##
##################################################################

targ_cmle <- all_pairs_nopc %>% 
  filter(cond_y == "CP-MLE") %>% 
  filter(cond_x != "CP-MLE")
targ_lpp <- all_pairs_nopc %>% 
  filter(cond_y == "PP-LT-C") %>% 
  filter(cond_x != "PP-LT-C", 
         cond_x != "CP-MLE")

targ_both <- bind_rows(targ_cmle, targ_lpp) %>% 
  filter(cond_x %in% sel_methods) %>% 
  droplevels 

## set
CUT_WIDTH_RHO <- 0.025
CUT_WIDTH_SE <- 0.02


### empirical plot
ppred_emp <- make_emp_biv_plot(targ_both, 
                               rhos_max, "Correlations (max)",
                               se_c, "SE (combined)") +
  facet_grid(~cond_label+cond_iv_label) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.1)) +
  coord_fixed(ratio = 1, xlim = c(0, 0.27), ylim = c(0, 0.95)) +
  theme(legend.position = "right") +
  scale_fill_gradientn(limits = c(0, 0.3), 
                         colours = c("darkgreen", "yellow", 
                                     "orange", "red", "darkred"), 
                         values = 
                           scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2, 0.3)),
                         name = "Absol.\ndev.", 
                         na.value = "transparent")

ggsave("figures_man/pred_empirical.pdf", plot = ppred_emp,
       width = 19, height = 13, units = "cm")

targ_both %>% 
  filter(se_c < CUT_WIDTH_SE, rhos_max < CUT_WIDTH_RHO) %>% 
  group_by(cond_label, cond_iv_label) %>% 
  summarise(
    m_abs_dev = mean(abs_dev),
    sd_abs_dev = sd(abs_dev),
    n = n()
  )

options(pillar.sigfig = 4)
targ_both %>% 
  filter(se_c < CUT_WIDTH_SE, rhos_max < CUT_WIDTH_RHO) %>% 
  group_by(cond_label, cond_iv_label) %>% 
  arrange(rhos_max) %>% 
  select(cond_label, cond_iv_label, abs_dev, rhos_max, se_c, 
         model2, dataset, parameter)
options(pillar.sigfig = 2)

# ## 3 smallest correlations and absolute deviations for SE < .02 and each pair.
# targ_both %>% 
#   filter(se_c < .02) %>% 
#   group_by(cond_label, cond_iv_label) %>% 
#   arrange(rhos_max) %>% 
#   select(cond_label, cond_iv_label, abs_dev, rhos_max, se_c) %>% 
#   slice(1:3)


ppred_gam <- make_gam_biv_plot(targ_both, 
                               rhos_max, "Correlations (max)",
                               se_c, "SE (combined)", 
                               return = "both") 

ppred_gam_plot <- ppred_gam$plot +
  facet_grid(~cond_label+cond_iv_label) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.1)) +
  coord_fixed(ratio = 1, xlim = c(0, 0.27), ylim = c(0, 0.95)) +
  theme(legend.position = "right") +
  scale_fill_gradientn(limits = c(0, 0.3), 
                         colours = c("darkgreen", "yellow", 
                                     "orange", "red", "darkred"), 
                         values = 
                           scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2, 0.3)),
                         name = "Mean\nabsol.\ndev.", 
                         na.value = "transparent")
## the warning is fine, within each panel the horizontal intervals are even

ggsave("figures_man/bigam.pdf",
       width = 19, height = 13, units = "cm")

ppred_legend <- get_legend(ppred_emp)

pl_left <- plot_grid(
  ppred_emp + theme(legend.position = "none", 
                    axis.text.x = element_blank(), 
                    axis.title.x = element_blank()),
  ppred_gam_plot + theme(legend.position = "none", 
                         strip.text = element_blank()),
  nrow = 2, rel_heights = c(1, 0.975)
  )

plot_out <- plot_grid(
  pl_left, ppred_legend, nrow = 1, rel_widths = c(1, 0.1)
)
ggsave("figures_man/ppred_both.pdf", width = 19, height = 24, units = "cm")

nd <- data.frame(se_c = 0, rhos_max = 0)

res_pred <- bind_cols(ppred_gam$gams, 
          map_dfr(ppred_gam$gams$res, 
                  ~as.data.frame(insight::get_predicted(., data = nd))))
res_pred %>% 
  select(-c(data:nd))


##################################################################
##        Selected Pairs of Target and Predicting Method        ##
##################################################################

rel_vars <- c("se_c", "y", "rel_weight", "rel_n_w", "rhos_max")
other_vars <- c("sd_emp_inv", "fungi_max", "p_hetero", "p_fit_x")
all_vars <- c(rel_vars, other_vars)
form_main <- paste0("abs_dev ~ ", paste(all_vars, collapse = " + "))

all_pairs_nopc_red <- all_pairs_nopc %>% 
  mutate(filtervar = create_filter_from_formula(form_main)) %>% 
  filter(filtervar)

1 - nrow(all_pairs_nopc_red)/nrow(all_pairs_nopc)

## see corresponding Report in docs folder for full overview
str(all_pairs)

targ_cmle_2 <- all_pairs_nopc_red %>% 
  filter(cond_y == "CP-MLE") %>% 
  filter(cond_x != "CP-MLE")
targ_lpp_2 <- all_pairs_nopc_red %>% 
  filter(cond_y == "PP-LT-C") %>% 
  filter(cond_x != "PP-LT-C", 
         cond_x != "CP-MLE")

targ_both_2 <- bind_rows(targ_cmle_2, targ_lpp_2) %>% 
  filter(cond_x %in% sel_methods) %>% 
  droplevels
ylab <- "Abs. deviation"

######


##----------------------------------------------------------------
##                  Absolute Mean Deviation Plot                 -
##----------------------------------------------------------------

targ_both_2 %>% 
  ggplot(aes(x = cond_x, y = abs_dev)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_violin(fill = "transparent", width = 1.1) +
  stat_summary(fun = mean, fun.max = mean, fun.min = mean, fatten = 0.9) +
  facet_grid(cols = vars(cond_label), scales = "free_x", space = "free_x") +
  labs(x = "Comparison method", y = ylab)
ggsave("figures_man/mad.pdf", 
       width = 16, height = 7, units = "cm", 
       dpi = 500)


##----------------------------------------------------------------
##                          RMSE table                           -
##----------------------------------------------------------------

rmse_tab_1 <- targ_both_2 %>% 
  group_by(cond_label, cond_iv_label) %>% 
  summarise(mean = mean(abs_dev), 
            sd = sd(abs_dev), 
            rmse = sqrt(mean((mean(abs_dev) - abs_dev)^2)), 
            sigma = sigma(lm(abs_dev ~ 1))) %>% 
  mutate(across(where(is.double), ~formatC(.x, digits = 3, format = "f")))
rmse_tab_1

rmse_tab_2 <- targ_both_2 %>% 
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

rmse_tab_3 <- targ_both_2 %>% 
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


##----------------------------------------------------------------
##                    Univariate Relationships                   -
##----------------------------------------------------------------


### parameter-level covariates
pest <- compare_continuous_covariate(data = targ_both_2, covariate = poly(y, 2), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Value of parameter estimate (reference method, quadratic)")

psec <- compare_continuous_covariate(data = targ_both_2, covariate = se_c, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (combined)")
psd <- compare_continuous_covariate(data = targ_both_2, covariate = sd_emp_inv, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Individual variability (SD)")
prho <- compare_continuous_covariate(data = targ_both_2, covariate = rhos_max, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter correlations (max)")
pfungi <- compare_continuous_covariate(data = targ_both_2, covariate = fungi_max, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter trade-offs (max)")
prelw <- compare_continuous_covariate(data = targ_both_2, covariate = log(rel_weight), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Relative information (log)")
preln <- compare_continuous_covariate(data = targ_both_2, covariate = log(rel_n_w), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Relative N (log)")

## data-set-level covariates
phetero <- compare_continuous_covariate(data = targ_both_2, covariate = log1p(p_hetero), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Hetereogeneity (log p + 1)")
pfit <- compare_continuous_covariate(data = targ_both_2, covariate = log1p(p_fit_x), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Model fit (comparison method, log p + 1)")

#breaks1 <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
#labels1 <- c("0", "0.2", "0.4", "0.6", "0.8", "1")
breaks1 <- c(0, 0.25, 0.5, 0.75, 1)
labels1 <- c("0", "0.25", "0.5", "0.75", "1")

library("cowplot")
puniv_good <- cowplot::plot_grid(
  psec + theme(axis.title.y = element_text(colour = "transparent")) +
    scale_x_continuous(breaks = c(0, 0.1, 0.2)),
  prho + theme(strip.text = element_blank(), 
               axis.title.y = element_text(colour = "transparent")), 
  pest + theme(strip.text = element_blank(), 
               axis.title.y = element_text(colour = "transparent")) + 
    scale_x_continuous(breaks = breaks1, labels = labels1),
  psd + theme(strip.text = element_blank(), 
              axis.title.y = element_text(colour = "transparent")),
  #pfungi + theme(strip.text = element_blank()), 
  ncol = 1, rel_heights = c(1.33, rep(1, 4))) +
  draw_text("Abs. deviation", angle = 90, x = 0.015)

ggsave("figures_man/univariate_good.png", plot = puniv_good, 
       width = 24.5, height = 21, units = "cm", 
       dpi = 500)

puniv_notgood <- cowplot::plot_grid(
   prelw + theme(axis.title.y = element_text(colour = "transparent")) +
   scale_x_continuous(breaks = breaks1, labels = labels1),  
  preln + theme(strip.text = element_blank(),
                axis.title.y = element_text(colour = "transparent")) +
     scale_x_continuous(breaks = c(0, 5000, 10000, 15000), 
                       labels = c("0", "5k", "10k", "15k")), 
  pfungi + theme(strip.text = element_blank()),
  phetero + theme(strip.text = element_blank(),
                  axis.title.y = element_text(colour = "transparent")) + 
    scale_x_continuous(breaks = breaks1, labels = labels1), 
  pfit + theme(strip.text = element_blank(),
               axis.title.y = element_text(colour = "transparent")) + 
    scale_x_continuous(breaks = breaks1, labels = labels1), 
  ncol = 1, rel_heights = c(1.35, rep(1, 4)))
ggsave("figures_man/univariate_notgood.png", plot = puniv_notgood, 
       width = 24.5, height = 25, units = "cm", 
       dpi = 500)


##----------------------------------------------------------------
##                      Multivariate Results                     -
##----------------------------------------------------------------

### See corresponding RMarkdown document:
### multivariate_relationships_with_abs_dev_rmse_noPC_2


##################################################################
##                          Discussion                          ##
##################################################################


all_pairs %>% 
  distinct(model, model2, dataset, parameter, cond_co, 
           parameter_o, condition, orig_condition, .keep_all = TRUE) %>% 
  filter(cond_x != cond_y) %>% 
  group_by(model, cond_x, cond_y) %>% 
  summarise(max = max(abs_dev)) %>% 
  arrange(max)
  

all_pairs %>% 
  distinct(model, model2, dataset, parameter, cond_co, 
           parameter_o, condition, orig_condition, .keep_all = TRUE) %>% 
  filter(cond_x != cond_y) %>% 
  filter(cond_x == "PP-LT-C") %>% 
  group_by(model, cond_x, cond_y) %>% 
  summarise(max = max(abs_dev)) %>% 
  arrange(max)

all_pairs %>% 
  distinct(model, model2, dataset, parameter, cond_co, 
           parameter_o, condition, orig_condition, .keep_all = TRUE) %>% 
  filter(cond_x != cond_y) %>% 
  #filter(cond_x == "PP-LT-C") %>% 
  group_by(model) %>% 
  arrange(desc(abs_dev)) %>% 
  slice(1:3) %>% 
  select(model:abs_dev, cond_x, cond_y, -dataset) %>% 
  print(n = Inf)



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


