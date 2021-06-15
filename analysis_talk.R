
library("afex")
library("DescTools")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))
source("fun_analysis_paper.R")
source("docs/fun_eai.R")

library("mgcv")
library("gratia")
INCLUDE_GAM <- TRUE

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
  mutate(prop_ns = rowMeans(cbind(prop_ns_nopb, prop_ns_nonpb, prop_ns_noasy, prop_ns_trait), na.rm = TRUE)) %>% 
  mutate(se_c = (se_x_w + se_y_w)/2) %>% 
  mutate(SAI = factor(if_else(model == "pc", "Not SAI", "SAI"), 
                      levels = c("SAI", "Not SAI"))) %>% 
  mutate(z_abs_dev = abs_dev/se_c) %>% 
  filter(parameter != "rm:g")
all_pairs <- all_pairs %>%
  mutate(cond_co = apply(cbind(as.character(cond_x), 
            as.character(cond_y)), 1, 
      function(x) paste(sort(x), collapse = " - ")))


##################################################################
##        Selected Pairs of Target and Predicting Method        ##
##################################################################


theme_set(theme_bw(base_size = 13) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))



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
sel_methods <- c("CP-MLE", "NP-MLE", "PP-B", "PP-LT-C")

targ_cmle <- all_pairs_red %>% 
  filter(cond_y == "CP-MLE") %>% 
  filter(cond_x != "CP-MLE")
targ_lpp <- all_pairs_red %>% 
  filter(cond_y == "PP-LT-C") %>% 
  filter(cond_x != "PP-LT-C")

targ_both <- bind_rows(targ_cmle, targ_lpp) %>% 
  filter(cond_x %in% sel_methods) %>% 
  droplevels %>% 
  mutate(cond_label = factor(paste("R:", cond_y), 
                             levels = c("R: CP-MLE", "R: PP-LT-C"))) %>% 
  mutate(cond_iv_label = factor(paste("C:", cond_x), levels = 
                                  paste(paste("C:", levels(cond_x)))))

ylab <- "Abs. deviation"


targ_both_nopc <- targ_both %>% 
  filter(model != "pc")



##----------------------------------------------------------------
##                      SAI deviation plot                       -
##----------------------------------------------------------------

all_pairs %>% 
  filter(cond_x %in% sel_methods, cond_y %in% sel_methods) %>% 
  filter(cond_y != "PP-B", 
         !(cond_x == "NP-MLE" & cond_y == "PP-LT-C"), 
         !(cond_x == "CP-MLE" & cond_y == "PP-LT-C"),
         !(cond_x == "CP-MLE" & cond_y == "NP-MLE")) %>% 
  ggplot(aes(x = SAI, y = abs_dev)) +
  #geom_boxplot(outlier.alpha = 0.1) +
  geom_violin() +
  stat_summary(color = "red") +
  facet_grid(cond_x2~ cond_y2, switch = "y", as.table = TRUE) +
  labs(y = ylab, x = "Structural Aggregation Invariance (SAI)")

ggsave("figures_talk/saiplot.png", 
       width = 15, height = 13, units = "cm", 
       dpi = 500)


##----------------------------------------------------------------
##                          Pairs Plot                           -
##----------------------------------------------------------------

plot_text <- all_pairs %>%
  filter(cond_x %in% sel_methods, cond_y %in% sel_methods) %>% 
  filter(model != "pc") %>% 
  group_by(cond_x2, cond_y2) %>% 
  summarise(ccc = substr(formatC(
    CCC(x, y, na.rm = TRUE)$rho.c$est, 
    digits = 2, format = "f"), 2, 4))
# plot_text$ccc

plot_rmse <- all_pairs %>%
  filter(cond_x %in% sel_methods, cond_y %in% sel_methods) %>% 
  filter(model != "pc") %>% 
  group_by(cond_x2, cond_y2) %>% 
  summarise(rmse = substr(formatC(
    sqrt(mean((x - y)^2)), 
    digits = 3, format = "f"), 2, 5))

ppairs <- all_pairs %>% 
  filter(cond_x %in% sel_methods, cond_y %in% sel_methods) %>% 
  filter(model != "pc") %>% 
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
ggsave("figures_talk/pairsplot_noPC.png", plot = ppairs,
       width = 17, height = 16, units = "cm", 
       dpi = 500)


##---------------------------------------------------------------
##                        EAI GAM Plots                         -
##---------------------------------------------------------------

make_gam_biv_plot(targ_both_nopc, 
                  rhos_max, "rel. correlations (max)",
                  se_c, "SE (combined)") +
   facet_grid(~cond_label+cond_iv_label)

ggsave("figures_talk/bigam.png",
       width = 21, height = 12, units = "cm", 
       dpi = 500)

##----------------------------------------------------------------
##        Zoom in: Effect of individual differences on EAI       -
##----------------------------------------------------------------

MIN_REL_SE <- 0.025
MIN_REL_RHO <- 0.15

targ_both_nopc %>% 
  filter(rhos_max < MIN_REL_RHO, se_c < MIN_REL_SE) %>%
  make_cor_plot2(x = prop_ns, y = abs_dev, 
                 filter = prop_ns < 0.4) +
  labs(x = "Proportion individual misfits (NP and PP-LT p < .05)", 
       y = ylab)

ggsave("figures_talk/corr.png",
       width = 18, height = 9, units = "cm", 
       dpi = 500)

##----------------------------------------------------------------
##                    Univariate Relationships                   -
##----------------------------------------------------------------

theme_set(theme_bw(base_size = 13) + 
            theme(legend.position="bottom"))

### parameter-level covariates
pest <- compare_continuous_covariate(data = targ_both_nopc, covariate = poly(y, 2), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Value of parameter estimate (reference method, quadratic)")
psex <- compare_continuous_covariate(data = targ_both_nopc, covariate = se_x_w, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (comparison method)")
psey <- compare_continuous_covariate(data = targ_both_nopc, covariate = se_y_w, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (reference method)")
psec <- compare_continuous_covariate(data = targ_both_nopc, covariate = se_c, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("SE (combined)")


psd <- compare_continuous_covariate(data = targ_both_nopc, covariate = sd_emp_inv, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Individual variability (SD)")
prho <- compare_continuous_covariate(data = targ_both_nopc, covariate = rhos_max, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter correlations (max)")
pfungi <- compare_continuous_covariate(data = targ_both_nopc, covariate = fungis_max, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Parameter trade-offs (max)")
prelw <- compare_continuous_covariate(data = targ_both_nopc, covariate = log(rel_weight), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Relative information (log)")
preln <- compare_continuous_covariate(data = targ_both_nopc, covariate = log(rel_n_w), 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Relative N (log)")

## data-set-level covariates
phetero <- compare_continuous_covariate(data = targ_both_nopc, 
                                        covariate = hetero_cohenw, 
                                     cond_label, cond_iv_label, ylab = ylab) +
  xlab("Hetereogeneity")

phetero_b <- targ_both_nopc %>% 
  mutate(chisq_hetero = if_else(chisq_hetero > 7500, 7500, chisq_hetero)) %>% 
  compare_continuous_covariate(covariate = chisq_hetero, 
                               cond_label, cond_iv_label, ylab = ylab) +
  xlab("Hetereogeneity")

pfit <- compare_continuous_covariate(data = targ_both_nopc, 
                                     covariate = fit_cohenw_x, 
                                     cond_label, cond_iv_label, ylab = ylab, 
                                     scales = "free_x") +
  xlab("Model fit (comparison method)")

pfit_b <- 
  compare_continuous_covariate(data = targ_both_nopc, covariate = stat_fit_x, 
                                     cond_label, cond_iv_label, ylab = ylab, 
                               scales = "free_x") +
  xlab("Model fit")

pifit <- 
  compare_continuous_covariate(data = targ_both_nopc, covariate = prop_ns, 
                                     cond_label, cond_iv_label, ylab = ylab, 
                               filter = prop_ns < 0.4) +
  xlab("Proportion individual misfits (NP and PP-LT p < .05)")


newbreaks <- function (n = 5, ...) 
{
    n_default <- n
    function(x, n = n_default) {
        x <- x[is.finite(x)]
        if (length(x) == 0) {
            return(numeric())
        }
        rng <- range(x)
        rng[2] <- rng[2]*0.8
        labeling::extended(rng[1], rng[2], n, ...)
    }
}

puniv_good <- cowplot::plot_grid(
  psec + scale_x_continuous(breaks = newbreaks(3)), 
  prho + theme(strip.text = element_blank()),
  pfungi + theme(strip.text = element_blank()),
  pifit + theme(strip.text = element_blank()) + 
    scale_x_continuous(breaks = newbreaks(3)), 
  pest + theme(strip.text = element_blank()) +
    scale_x_continuous(breaks = newbreaks(3)), 
  ncol = 1, rel_heights = c(1.3, rep(1, 4)))
ggsave("figures_talk/univariate_good_noPC.png", plot = puniv_good, 
       width = 28, height = 25, units = "cm", 
       dpi = 500)



puniv_notgood <- cowplot::plot_grid(
  prelw + scale_x_continuous(breaks = newbreaks(3)), 
  preln + theme(strip.text = element_blank()), 
  phetero + theme(strip.text = element_blank()), 
  pfit + theme(strip.text = element_blank()) + 
    scale_x_continuous(
      breaks = newbreaks(3)), 
  psd + theme(strip.text = element_blank()),
  ncol = 1, rel_heights = c(1.3, rep(1, 4)))
ggsave("figures_talk/univariate_notgood_noPC.png", plot = puniv_notgood, 
       width = 28, height = 25, units = "cm", 
       dpi = 500)

