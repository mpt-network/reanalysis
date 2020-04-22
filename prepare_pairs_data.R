
library("MPTmultiverse")
library("tidyverse")
library("tidylog")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom"))
tmpe <- new.env()
source("fun_prep.R", local = tmpe)
check_set <- get(x = "check_set",envir = tmpe)
library("DescTools") # for CCC

load("combined_results.RData")
source("fun_analysis.R")
covariates <- read_csv("covariates.csv")

lev_mod <- levels(dm$model)
lev_mod2 <- levels(dm$model2)

dm <- left_join(dm, covariates) %>% 
  mutate(model = factor(model, levels = lev_mod),
         model2 = factor(model2, levels = lev_mod2))

#### condition names

dm %>% 
  unnest(est_group) %>% 
  group_by(model, model2, dataset, orig_condition) %>% 
  summarise() %>% 
  write_csv("conditon_names.csv")

#####

dm$use <- NA
for (i in seq_len(nrow(dm))) {
  dm[i, "use"] <- check_fit(dm[i,])
}
prop.table(table(dm$use, useNA = "ifany"))

levels(dm$inter)
levels(dm$model)

exclude_methods <- c("LC PP") ## could be a value of inter
exclude_models <- c("")  

results <- dm %>% 
  filter(use) %>% 
  filter(!(inter %in% exclude_methods)) %>% 
  filter(!(model %in% exclude_models)) %>% 
  droplevels

##----------------------------------------------------------------
##                          Covariates                           -
##----------------------------------------------------------------

results$npar <- NA_integer_

for (i in seq_len(nrow(results))) {
  results$npar[[i]] <- results$est_group[[i]] %>% 
    filter(condition == condition[1]) %>% 
    nrow()
}

results$npar_c <- results$npar - 6

results %>% 
  select(npar, npar_c) %>% 
  psych::describe()

# results$est_rho[1:10]
# 
# results$fungibility[1:10]


##---------------------------------------------------------------
##                        Failure Rates                         -
##---------------------------------------------------------------


## get_convergence(results, model, split = TRUE) ##fails

convergence <- get_convergence(results, split = FALSE)
convergence_data <- results %>%
    group_split(dataset) %>% 
    {suppressWarnings(map_dfr(., check_set))}
convergence_data$sum <- rowSums(convergence_data[,-(1:2)])

table(convergence_data$sum) %>% 
  prop.table
#          7          8          9 
# 0.03012048 0.10843373 0.86144578 

results %>% 
  #filter(model != "htsm") %>% 
  get_convergence(model, split = FALSE)
# # A tibble: 11 x 2
#    key           value
#    <chr>         <dbl>
#  1 Comp_MR_MLE   1    
#  2 Comp_TB_ss    1    
#  3 No_MR_MLE     1    
#  4 No_MR_NPB     0.982
#  5 No_MR_PB      1    
#  6 No_TB_ss      0.964
#  7 PP_LC_lc      0    
#  8 PP_TB_beta    0.988
#  9 PP_TB_beta++  0    
# 10 PP_TB_trait   0.946
# 11 PP_TB_trait_u 0.952

convergence %>% 
  filter(value != 0) %>% 
  mutate(method = factor(key, levels = c(
    "Comp_MR_MLE", "No_MR_MLE", "No_MR_NPB", "No_MR_PB", 
    "Comp_TB_ss", "No_TB_ss", "PP_TB_beta", "PP_TB_trait_u", "PP_TB_trait"), 
    labels = c("Comp MLE", "No asy", "No PB", "No NPB", 
               "Comp Bayes", "No Bayes",  
               "Beta PP", "Trait_u PP", "Trait PP"))) %>% 
  mutate(failure = 1 - value) %>% 
  ggplot(aes(x = failure, y = method)) +
  geom_point(size = 5) +
  labs(y = "Method", x = "Failure Rate") +
  coord_cartesian(xlim = c(0, 0.055)) +
  theme_bw(base_size = 25) + 
  scale_x_continuous(breaks = seq(0, 0.05, by = 0.01), 
                     labels = scales::label_percent(accuracy = 1))
ggsave("figures/failure.png", width = 16, height = 18, units = "cm",
       dpi = 500)


##################################################################
##                         CCC Analysis                         ##
##################################################################

##----------------------------------------------------------------
##            Combine individual parameter estimates             -
##----------------------------------------------------------------

all_pars <- results %>% 
  select(-orig_condition) %>% 
  unnest(est_group) %>% 
  mutate(parameter_only = parameter) %>% 
  mutate(parameter = factor(paste0(model2, ":", parameter))) %>% 
  mutate(parameter_o = factor(paste0(model2, ":", orig_parameter))) %>% 
  droplevels()

##---------------------------------------------------------------
##                        More Covariates                       -
##---------------------------------------------------------------

## add covariates:

# all_pars <- all_pars %>% 
#   left_join(covariates_par) %>% 
#   left_join(covariates_nopar)


any(is.na(all_pars$core))

length(unique(all_pars$parameter))

options("tidylog.display" = list())
all_pars <- all_pars %>% 
  group_by(model, dataset, inter, orig_model, orig_condition) %>% 
  mutate(rel_par_weight = 
           get_rel_par_weight(
             parameter = orig_parameter, 
             est = est, 
             model_exp = model_exp, rel_tree = rel_tree, 
             orig_condition = orig_condition), 
         rel_n = get_weight_n(
             parameter = orig_parameter, 
             est = est, 
             model_exp = model_exp, data_tree = data_tree, 
             orig_condition = orig_condition)) %>% 
  ungroup
options("tidylog.display" = NULL)

# all_pars %>% 
#   filter(is.na(rel_par_weight)) %>% 
#   select(model, model2, dataset, rel_par_weight) %>% 
#   print(n = Inf)

### missing relative weights
all_pars %>% 
  filter(is.na(rel_par_weight)) %>% 
  group_by(model, dataset) %>% 
  tally %>% 
  print(n = Inf)

all_pars %>% 
  filter(is.na(rel_n)) %>% 
  group_by(model, dataset, orig_parameter) %>% 
  tally 


### N-participants (not used, we focus on SE instead)
# uc <- unique(all_pars$condition)
# uoc <- unique(n_participants$orig_condition)
# uc[ !(uc %in% uoc)]
# 
# n_participants2 <- all_pars %>% 
#   select(model, dataset, condition, model2, orig_condition) %>% 
#   left_join(n_participants)

core_pars <- all_pars %>% 
  filter(core) %>% 
  droplevels


#### covariates only

length(unique(core_pars$parameter))


fungibility <- results %>% 
  filter(inter == "Trait PP") %>% 
  unnest(fungibility) %>% 
  gather("which_par", "parameter", parameter1, parameter2) %>% 
  mutate(parameter_o = factor(paste0(model2, ":", parameter))) %>%
  group_by(model, dataset, inter, model2, 
           npar, npar_c, condition, parameter_o) %>% 
  summarise(fungi = max(abs(correlation), na.rm = TRUE)) %>% 
  ungroup() %>% 
  rename(orig_condition = condition)

correlation <- results %>% 
  filter(inter == "Trait PP") %>% 
  unnest(est_rho) %>% 
  gather("which_par", "parameter", parameter1, parameter2) %>% 
  group_by(model, dataset, inter, model2, 
           npar, npar_c, condition, parameter) %>% 
  summarise(rho = mean(abs(est))) %>% 
  ungroup() %>% 
  rename(orig_condition = condition) %>% 
  mutate(parameter_o = factor(paste0(model2, ":", parameter))) %>% 
  select(-parameter)

## covariates on parameters basis
covariates_par <- left_join(fungibility, correlation) 
## Note: models where Trait PP failed, do not have values here


p_vals <- results %>% 
  unnest(gof_group) %>% 
  filter(focus == "mean") %>% 
  select(model, dataset, inter, model2, condition, p, npar) %>% 
  rename(orig_condition = condition) %>% 
  rename(p_fit = p)

hetero_np <- results %>% 
  filter(inter == "Comp MLE") %>% 
  unnest(test_homogeneity) %>% 
  select(model, dataset, model2, condition:p) %>% 
  rename(orig_condition = condition) %>% 
  rename(p_hetero = p)

## covariates not relative to parameter
covariates_nopar <- left_join(p_vals, hetero_np)

## should be of length 0
covariates_nopar %>% 
  filter(is.na(p_fit) | is.na(p_hetero)) %>% 
  print(n = Inf)

### all data in covariates nopar:
unique(select(all_pars, model, dataset, inter, orig_condition)) %>% 
  anti_join(unique(select(covariates_nopar, model, dataset, inter, orig_condition)))

### should be of length 0!

### add parameters in which PP is missing
missing_par_covariates <- all_pars %>% 
  filter(inter == "Comp MLE") %>%
  select(model, model2, dataset, parameter_o, npar, npar_c, orig_condition) %>% 
  anti_join(select(covariates_par, -inter))

all_par_covariates <- full_join(select(covariates_par, -inter), 
                                missing_par_covariates)

all_par_covariates %>% 
  filter(is.na(orig_condition))

all_par_covariates %>% 
  filter(is.na(model2))

all_par_covariates %>% 
  group_by(model2) %>% 
  summarise(mean(is.na(fungi)))

all_par_covariates %>% 
  filter(model2 == "pd_e") %>% 
  print(n = Inf) ## still contains non core parameters


xx1 <- all_par_covariates %>% 
  select(model, dataset, orig_condition, parameter_o) %>% 
  unique 
xx2 <- unique(select(all_pars, model, dataset, orig_condition, parameter_o))
anti_join(xx2, xx1)
### should be of length 0!

covariates <- full_join(select(all_par_covariates, -npar), 
                        covariates_nopar)

### add real parameter names
covariates <- covariates %>% 
  left_join(select(all_pars, model, dataset, orig_condition, inter, parameter_o, 
                   parameter, condition, rel_par_weight, rel_n, se))

##----------------------------------------------------------------
##          Calculate CCC of all pairwise combinations           -
##----------------------------------------------------------------

#pairs <- combn(sort(levels(core_pars$inter)), 2)
pairs <- combn(levels(core_pars$inter), 2)

all_pars_l <- vector("list", ncol(pairs))

for (i in seq_len(ncol(pairs))) {
  tmp_dat <- core_pars %>% 
    filter(inter %in% pairs[,i]) %>% 
    mutate(inter = factor(inter, 
                           levels = unique(as.character(unlist(pairs[,i]))))) %>% 
    select(model, model2, dataset, condition, parameter, est, inter) %>% 
    spread(key = inter, value = est)
  colnames(tmp_dat)[(ncol(tmp_dat)-1):ncol(tmp_dat)] <- c("x", "y")
  tmp_dat$cond_x <- pairs[1,i]
  tmp_dat$cond_y <- pairs[2,i]
  all_pars_l[[i]] <- tmp_dat 
}
all_pars2 <- bind_rows(all_pars_l)
#all_pars2 <- left_join(all_pars2, trials)
all_pars2 <- all_pars2 %>% 
  mutate(
    cond_x = factor(cond_x, levels = 
                      levels(all_pars$inter)[unique(all_pars2$cond_x) %in% levels(all_pars$inter)]),
    cond_y = factor(cond_y, levels = 
                      levels(all_pars$inter)[unique(all_pars2$cond_y) %in% levels(all_pars$inter)])
    )

all_pars3 <- all_pars2 %>%
  # filter(restriction == "rrest", parameter %in% c("Do", "Dn", "g")) %>%
  # mutate(parameter = factor(parameter,
  #                           levels = c("Do", "Dn", "g"))) %>%
  mutate(cond_x2 = cond_y,
         cond_y2 = cond_x)



# all_ccc <- all_pars2 %>%
#   group_by(cond_x, cond_y) %>% 
#   summarise(ccc = CCC(x, y, na.rm = TRUE)$rho.c$est) %>% 
#   ungroup() %>% 
#   mutate(id = seq_len(nrow(.)))
# 
# ## plot of CCC for "rrest" variant across parameters
# all_ccc %>% 
# ggplot(aes(x = cond_x, y = ccc)) +
#   ggbeeswarm::geom_beeswarm(aes(color = cond_y), 
#                             cex = 0.6, alpha = 0.9, #color = "darkgrey", 
#                             size = 1.5) +
#   stat_summary(fun.data = mean_se) 
# ggsave("cccs1.png", width = 28, height = 18, units = "cm", 
#        dpi = 500)

##----------------------------------------------------------------
##                      Selected CCC Plots                       -
##----------------------------------------------------------------

# all_pars_a3 <- all_pars_a3 %>% 
#   filter(!(model %in% "pm"))

# dput(levels(all_pars_a3$cond_x))
# dput(levels(all_pars_a3$cond_y))

all_pars3 %>% 
   filter(cond_x == "Comp MLE", cond_y == "Comp Bayes")


p1 <- plot_pair(all_pars3, "Comp MLE", "No PB", FALSE) +
  labs(x = "Complete Pooling MLE", y = "No Pooling: P-Bootstrap")

p2 <- plot_pair(all_pars3, "Comp MLE", "Trait PP", FALSE) +
  labs(x = "Complete Pooling MLE", y = "Bayesian Latent Trait")

p3 <- plot_pair(all_pars3, "No PB", "Trait PP", FALSE) +
  labs(x = "No Pooling: P-Bootstrap", y = "Bayesian Latent Trait")

cowplot::plot_grid(p2, p1, p3, nrow = 1)
ggsave("figures/pairsplot_selected_1.png", width = 36, height = 12, 
       units = "cm", dpi = 500)

cowplot::plot_grid(
  p2 + geom_smooth(), 
  p1 + geom_smooth(), 
  p3 + geom_smooth(), 
  nrow = 1)
ggsave("figures/pairsplot_selected_2.png", width = 36, height = 12, 
       units = "cm", dpi = 500)

plot_pair(all_pars3, "No asy", "Trait PP", FALSE) +
  geom_smooth()

all_pars3 %>% 
  filter(cond_x == "Comp MLE", cond_y == "Trait PP") %>% 
  mutate(adev = abs(x - y)) %>% 
  arrange(desc(adev)) %>% 
  print(n = 40)

all_pars3 %>% 
  filter(cond_x == "Comp MLE", cond_y == "No PB") %>% 
  mutate(adev = abs(x - y)) %>% 
  arrange(desc(adev))

all_pars3 %>% 
  filter(cond_x == "No PB", cond_y == "Trait PP") %>% 
  mutate(adev = abs(x - y)) %>% 
  arrange(desc(adev)) %>% 
  print(n = 20)


### pairs plot without pair clustering model

p1 <- all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  plot_pair("Comp MLE", "No PB", FALSE) +
  labs(x = "Complete Pooling MLE", y = "No Pooling: P-Bootstrap")

p2 <- all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  plot_pair("Comp MLE", "Trait PP", FALSE) +
  labs(x = "Complete Pooling MLE", y = "Bayesian Latent Trait")

p3 <- all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  plot_pair("No PB", "Trait PP", FALSE) +
  labs(x = "No Pooling: P-Bootstrap", y = "Bayesian Latent Trait")

cowplot::plot_grid(p2, p1, p3, nrow = 1)
ggsave("figures/pairsplot_selected_1-no-PC.png", width = 36, height = 12, 
       units = "cm", dpi = 500)

all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  filter(cond_x == "Comp MLE", cond_y == "Trait PP") %>% 
  mutate(adev = abs(x - y)) %>% 
  arrange(desc(adev))

all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  filter(cond_x == "Comp MLE", cond_y == "No PB") %>% 
  mutate(adev = abs(x - y)) %>% 
  arrange(desc(adev))

all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  filter(cond_x == "No PB", cond_y == "Trait PP") %>% 
  mutate(adev = abs(x - y)) %>% 
  arrange(desc(adev)) %>% 
  print(n = 20)

all_pars3 %>% 
  filter(!(model == "pc")) %>% 
  droplevels() %>% 
  {length(unique(.$parameter))}

all_pars3 %>% 
  filter(!(model == "pc")) %>% 
   filter(cond_x == "Comp MLE", cond_y == "Comp Bayes")

##----------------------------------------------------------------
##                        Pairs Plot All                         -
##----------------------------------------------------------------


pars_plot <- all_pars3 

plot_text <- pars_plot %>%
  group_by(cond_x2, cond_y2) %>% 
  summarise(ccc = format(
    CCC(x, y, na.rm = TRUE)$rho.c$est, 
    digits = 2))

pars_plot %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.15) + #aes(size = trials)
  facet_grid(cond_x2~ cond_y2, switch = "both", as.table = FALSE) +
  # geom_text(data=plot_text,
  #           aes(x = 0.2, y = 0.9, label=ccc), 
  #           parse = TRUE, inherit.aes=FALSE, size = 5) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_size(range = c(0.5, 2.5)) +
  labs(x = "", y = "") +
  theme(legend.position = "none")
ggsave("figures/pairsplot_all.png", width = 33, height = 28, units = "cm", 
       dpi = 500)

pars_plot %>% 
  filter(cond_x == "Comp Bayes") %>% 
  filter(x < 0.05, y > 0.5) %>% 
  print(n = 1e5)

#################

for (i in levels(pars_plot$model)) {
  pars_plot <- pars_plot %>% 
    mutate(Selected = factor(if_else(model == i, toupper(i), "other"), 
                             levels = c(toupper(i), "other"))) %>% 
    arrange(desc(Selected))
  
  pars_plot %>% 
    ggplot(aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(color = Selected), alpha = 0.33) + #aes(size = trials)
    facet_grid(cond_x2~ cond_y2, switch = "both", as.table = FALSE) +
    # geom_text(data=plot_text,
    #           aes(x = 0.2, y = 0.9, label=ccc), 
    #           parse = TRUE, inherit.aes=FALSE, size = 5) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_size(range = c(0.5, 2.5)) +
    labs(x = "", y = "") +
    theme(legend.position = "none") +
    scale_colour_manual(values = c("red", "black"))
  ggsave(paste0("figures/pairsplot_all_", i, ".png"), 
         width = 33, height = 28, units = "cm", dpi = 500)
  
}

for (i in levels(pars_plot$model)) {
  pars_plot <- pars_plot %>% 
    mutate(Selected = factor(if_else(model == i, toupper(i), "other"), 
                             levels = c(toupper(i), "other"))) %>% 
    arrange(Selected)
  
  pars_plot %>% 
    ggplot(aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(color = Selected), alpha = 0.33) + #aes(size = trials)
    facet_grid(cond_x2~ cond_y2, switch = "both", as.table = FALSE) +
    # geom_text(data=plot_text,
    #           aes(x = 0.2, y = 0.9, label=ccc), 
    #           parse = TRUE, inherit.aes=FALSE, size = 5) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_size(range = c(0.5, 2.5)) +
    labs(x = "", y = "") +
    theme(legend.position = "none") +
    scale_colour_manual(values = c("red", "black"))
  ggsave(paste0("figures/pairsplot_all_", i, "_2.png"), 
         width = 33, height = 28, units = "cm", dpi = 500)
}


#################################################################
##                  All Pairwise Combinations                  ##
#################################################################

pairs <- expand.grid(cond_x = sort(levels(core_pars$inter)), 
                     cond_y = sort(levels(core_pars$inter)))



all_pars_l2 <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  tmp_dat <- core_pars %>% 
    filter(core) %>% 
    filter(inter %in% unlist(pairs[i,])) %>% 
    mutate(inter2 = factor(inter, 
                           levels = unique(as.character(unlist(pairs[i,]))))) %>% 
    select(model, model2, dataset, condition, parameter, est, inter2) %>% 
    spread(key = inter2, value = est)
  if ( (colnames(tmp_dat)[ncol(tmp_dat)-1]) == "parameter" ) {
    tmp_dat$y <- tmp_dat[[ncol(tmp_dat)]]
  }
  colnames(tmp_dat)[(ncol(tmp_dat)-1):ncol(tmp_dat)] <- c("x", "y")
  tmp_dat$cond_x <- pairs[i, 1]
  tmp_dat$cond_y <- pairs[i, 2]
  all_pars_l2[[i]] <- tmp_dat 
}
all_pars_a2 <- bind_rows(all_pars_l2)
#all_pars2 <- left_join(all_pars2, trials)
all_pars_a2 <- all_pars_a2 %>% 
  mutate(
    cond_x = factor(cond_x, levels = 
                      levels(results$inter)),
    cond_y = factor(cond_y, levels = 
                      levels(results$inter))
    )
all_pars_a3 <- all_pars_a2 %>%
  # filter(restriction == "rrest", parameter %in% c("Do", "Dn", "g")) %>%
  # mutate(parameter = factor(parameter,
  #                           levels = c("Do", "Dn", "g"))) %>%
  mutate(cond_x2 = cond_y,
         cond_y2 = cond_x) %>% 
  filter(!(str_detect(dataset, "Smith et al 2014")))

all_pars_a3 %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.15) + #aes(size = trials)
  facet_grid(cond_x2~ cond_y2, switch = "both", as.table = FALSE) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_size(range = c(0.5, 2.5)) +
  geom_smooth() +
  labs(x = "", y = "") +
  theme(legend.position = "none")
ggsave("figures/pairsplot_no-all.png", width = 33, height = 28, 
       units = "cm", dpi = 500)


##################################################################
##                      Absolute Deviation                      ##
##################################################################

all_pars_a4 <- all_pars_a3 %>% 
  mutate(adev = abs(x - y)) %>% 
  filter(!(cond_x == cond_y)) %>% 
  filter(!is.na(x), !is.na(y)) %>% 
  mutate(x1 = x - mean(x)) %>% 
  mutate(x2 = x1^2) %>% 
  droplevels()

unique(all_pars_a4$parameter)
unique(covariates$parameter)

all_pars_a4 <- covariates %>% 
  filter(parameter %in% levels(all_pars_a4$parameter)) %>% 
  droplevels %>% 
  rename(cond_x = inter) %>% 
  right_join(all_pars_a4) %>% 
  mutate(log1p_fit_x = log1p(p_fit),
         logp_fit_x = log(p_fit),
         log1p_hetero = log1p(p_hetero),
         logp_hetero = log(p_hetero)) %>% 
  rename(rel_par_weight_x = rel_par_weight,
         rel_n_x = rel_n, 
         p_fit_x = p_fit, 
         se_x = se) 

all_pars_a4 <- covariates %>% 
  filter(parameter %in% levels(all_pars_a4$parameter)) %>% 
  droplevels %>% 
  rename(cond_y = inter) %>% 
  select(c("model", "dataset", "model2", "cond_y", "parameter", "condition"), 
         rel_par_weight, rel_n, p_fit, se) %>% 
  right_join(all_pars_a4) %>% 
  mutate(log1p_fit_y = log1p(p_fit),
         logp_fit_y = log(p_fit)) %>% 
  rename(rel_par_weight_y = rel_par_weight,
         rel_n_y = rel_n, 
         p_fit_y = p_fit, 
         se_y = se)

# all_pars_a4 %>% 
#   group_by(model2) %>% 
#   summarise(mean(is.na(fungi)))
# 
# all_pars_a4 %>% 
#   summarise(mean(is.na(fungi)))

# all_pars_a4 %>% 
#   select(model, model2, dataset, adev, cond_x, x, cond_y, y, 
#          parameter, se_x, se_y, rho)

par_pairs <- all_pars_a4

save(par_pairs, file = "all_pairs_core.RData")
save(par_pairs, file = "all_pairs_core_OLD.RData", version = 2)

