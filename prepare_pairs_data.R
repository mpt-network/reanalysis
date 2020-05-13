
CURDIR <- ""
try(CURDIR <- dirname(rstudioapi::getActiveDocumentContext()$path))
if (!file.exists(CURDIR)) try(CURDIR <- dirname(parent.frame(2)$ofile))
print(CURDIR)
setwd(CURDIR)

library(checkpoint)
checkpoint("2020-04-25", R.version = "4.0.0")
## do not compile from source on windows!

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

## combine with other data
lev_mod <- levels(dm$model)
lev_mod2 <- levels(dm$model2)

#### condition names

# dm %>% 
#   unnest(est_group) %>% 
#   group_by(model, model2, dataset, orig_condition) %>% 
#   summarise() %>% 
#   write_csv("conditon_names.csv")

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
options("tidylog.display" = list())
results$npar <- NA_integer_

for (i in seq_len(nrow(results))) {
  results$npar[[i]] <- results$est_group[[i]] %>% 
    filter(condition == condition[1]) %>% 
    nrow()
}
options("tidylog.display" = NULL)

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


##----------------------------------------------------------------
##            Combine individual parameter estimates             -
##----------------------------------------------------------------

all_pars <- results %>% 
  unnest(est_group) %>% 
  mutate(parameter_only = parameter) %>% 
  mutate(parameter = factor(paste0(model2, ":", parameter))) %>% 
  mutate(parameter_o = factor(paste0(model2, ":", orig_parameter))) %>% 
  droplevels()

### check if covariates exists for all variables (and not more)
### should be TRUE
stopifnot(all(unique(covariates$dataset) %in% unique(all_pars$dataset)))
stopifnot(all(unique(all_pars$dataset) %in% unique(covariates$dataset)))

all_pars <- left_join(all_pars, covariates) %>%
  mutate(model = factor(model, levels = lev_mod),
         model2 = factor(model2, levels = lev_mod2))


## should be length 0
all_pars %>% 
  filter(is.na(population)) %>% 
  select(model, dataset, orig_condition) %>% 
  unique()

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
  group_by(model, dataset, inter) %>% 
  tally 

all_pars %>% 
  filter(is.na(est)) %>% 
  group_by(model, dataset, inter) %>% 
  tally()


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

hetero_empirical <- results %>% 
  filter(inter == "Trait PP") %>% 
  unnest(est_indiv) %>% 
  group_by(model, dataset, model2, orig_condition, 
           orig_parameter, parameter) %>% 
  summarise(sd_emp = sd(est)) %>%
  ungroup %>% 
  mutate(parameter_o = factor(paste0(model2, ":", orig_parameter))) %>% 
  select(-parameter, -orig_parameter)
mean(is.na(hetero_empirical$sd_emp))

## covariates on parameters basis
covariates_par <- left_join(fungibility, correlation) 
covariates_par <- left_join(covariates_par, hetero_empirical)

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

stopifnot(all(p_vals$orig_condition %in% hetero_np$orig_condition))
stopifnot(all(hetero_np$orig_condition %in% p_vals$orig_condition))

## covariates not relative to parameter
covariates_nopar <- left_join(p_vals, hetero_np)

## should be of length 0
covariates_nopar %>% 
  filter(is.na(p_fit) | is.na(p_hetero)) %>% 
  print(n = Inf)

### all data in covariates nopar:
unique(select(all_pars, model, dataset, inter, orig_condition)) %>% 
  anti_join(unique(select(covariates_nopar, model, dataset, inter, 
                          orig_condition)))
### should be of length 0!

### add parameters in which PP is missing
missing_par_covariates <- all_pars %>% 
  filter(inter == "Comp MLE") %>%
  select(model, model2, dataset, parameter_o, npar, npar_c, orig_condition) %>% 
  anti_join(select(covariates_par, -inter))

missing_par_covariates %>% 
  group_by(model, dataset) %>% 
  tally()

all_par_covariates <- full_join(select(covariates_par, -inter), 
                                missing_par_covariates)

all_par_covariates %>% 
  filter(is.na(orig_condition))

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

# tt1 <- select(all_par_covariates, -npar) %>% 
#   filter(dataset == "A2013")
# tt2 <- covariates_nopar %>% 
#   filter(dataset == "A2013")
# 
# covariates <- full_join(tt1, tt2)

covariates <- full_join(select(all_par_covariates, -npar), 
                        covariates_nopar)

# nrow(covariates_nopar)
# nrow(all_par_covariates)

### add real parameter names
covariates <- covariates %>% 
  left_join(select(all_pars, model, dataset, orig_condition, inter, parameter_o, 
                   parameter, condition, rel_par_weight, rel_n, se))

stopifnot(nrow(covariates) == nrow(all_pars))

xxx <- covariates %>% 
  anti_join(select(all_pars, model, dataset, orig_condition, inter, parameter_o, 
                   parameter, condition, rel_par_weight, rel_n, se)) 
## should be zero...
xxx %>% 
  select(dataset, parameter_o)


### add external covariates

## external covariate files
ext_covariates <- read_csv("covariates.csv", 
                       locale = readr::locale(encoding = "UTF-8")) %>% 
  mutate(population = factor(population), 
         sci_goal = factor(
           `parameter estimation?`, 
           levels = c("yes", "no"), 
           labels = c("estimation", "model_comparison")
         ),
         model = factor(model, levels = levels(covariates$model)),
         model2 = factor(model2, levels = levels(covariates$model2))
  ) %>% 
  select(-`parameter estimation?`)
# str(ext_covariates)

## check if all variables are filled 
### should be FALSE
any(is.na(ext_covariates$population))
any(is.na(ext_covariates$sci_goal))

covariates <- left_join(covariates, ext_covariates)

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


##################################################################
##                      Absolute Deviation                      ##
##################################################################

all_pairs <- all_pars_a2 %>% 
  mutate(abs_dev = abs(x - y)) %>% 
  filter(!(cond_x == cond_y)) %>% 
  filter(!is.na(x), !is.na(y)) %>% 
  #mutate(x_centered = x - mean(x)) %>% 
  #mutate(x_centered_pow2 = x_centered^2) %>% 
  droplevels()

any(is.na(all_pairs$abs_dev))

# unique(all_pars_a4$parameter)
# unique(covariates$parameter)

# covariates %>%
#   filter(parameter %in% levels(all_pairs$parameter)) %>%
#   droplevels %>%
#   rename(cond_x = inter) %>%
#   anti_join(all_pairs) %>% 
#   select(model:parameter_o, cond_x)
## parameter is NA for "No PB"

all_pairs <- covariates %>% 
  filter(parameter %in% levels(all_pairs$parameter)) %>% 
  droplevels %>% 
  rename(cond_x = inter) %>% 
  right_join(all_pairs) %>% 
  mutate(log1p_fit_x = log1p(p_fit),
         logp_fit_x = log(p_fit),
         log1p_hetero = log1p(p_hetero),
         logp_hetero = log(p_hetero)) %>% 
  rename(rel_par_weight_x = rel_par_weight,
         rel_n_x = rel_n, 
         p_fit_x = p_fit, 
         se_x = se) 
any(is.na(all_pairs$abs_dev))

all_pairs <- covariates %>% 
  filter(parameter %in% levels(all_pairs$parameter)) %>% 
  droplevels %>% 
  rename(cond_y = inter) %>% 
  select(c("model", "dataset", "model2", "cond_y", "parameter", "condition"), 
         rel_par_weight, rel_n, p_fit, se) %>% 
  right_join(all_pairs) %>% 
  mutate(log1p_fit_y = log1p(p_fit),
         logp_fit_y = log(p_fit)) %>% 
  rename(rel_par_weight_y = rel_par_weight,
         rel_n_y = rel_n, 
         p_fit_y = p_fit, 
         se_y = se)

all_pairs_full <- all_pairs %>% 
  mutate(dataset = factor(dataset):model2, 
         condition = factor(condition):model2) %>% 
  droplevels

str(all_pairs_full)

all_pairs <- all_pairs_full %>% 
  select(model, model2, dataset, parameter, 
         abs_dev, x, se_x, cond_x, y, se_y, cond_y, 
         logp_hetero, log1p_hetero, ## Heterogeneity (non-parameteric)
         sd_emp, ## Hetereogeneity across parameter based on average of 
                  ## partial-pooling (i.e., empirical SD of individual latent 
                  ## trait parameters.
         rho, ## Rho (average correlation for each parameter)
         fungi, ##  Fungibility/across-chain correlations 
         log1p_fit_x, logp_fit_x, log1p_fit_y, logp_fit_y, ## model fit
         rel_par_weight_x, rel_par_weight_y, ## Absolute parameter values 
         ## product of previous branches. Relative information available 
         npar, ## number of parameters
         population, sci_goal, ### external covariates
         condition, orig_condition, parameter_o
         )

str(all_pairs)

# all_pars_a4 %>% 
#   group_by(model2) %>% 
#   summarise(mean(is.na(fungi)))
# 
# all_pars_a4 %>% 
#   summarise(mean(is.na(fungi)))

# all_pars_a4 %>% 
#   select(model, model2, dataset, adev, cond_x, x, cond_y, y, 
#          parameter, se_x, se_y, rho)

save(all_pairs, file = "all_pairs_core.RData", compress = "xz")
save(all_pairs, file = "all_pairs_core_OLD.RData", compress = "xz", 
     version = 2)

save(all_pairs_full, file = "all_pairs_core_full.RData", compress = "xz")
save(all_pairs_full, file = "all_pairs_core_full_OLD.RData", compress = "xz", 
     version = 2)

