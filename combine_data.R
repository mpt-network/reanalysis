
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## https://ukgovdatascience.github.io/rap-website/article-dependency-and-reproducibility.html
library(checkpoint)
checkpoint("2020-04-25", R.version = "4.0.0")
## do not compile from source on windows!

library("MPTmultiverse")
library("tidyverse")
# library("tidylog")
source("fun_prep.R")

##################################################################
##                      Prospective Memory                      ##
##################################################################
## Nina Arnhold & Sebastian Horn

## latent class missing

f_pm <- list.files("data/prospective-memory/", full.names = TRUE)
fi_pm <- f_pm[str_detect(f_pm, "_info.RData$")]
f_pm <- f_pm[!str_detect(f_pm, "_info.RData$")]
pm <- load_combine(f_pm) %>% 
  mutate(dataset = str_remove(dataset, ".+\\/")) %>% 
  mutate(model = str_remove(model, ".+\\/")) %>% 
  mutate(orig_model = str_remove(orig_model, ".+\\/"))
get_info_df(pm)

check_beta(pm) 

pm$model <- "pm"
pm$model2 <- "pm"

pm <- check_double_fit(pm)

## change within-subject variable names:
### move addition to variable name to condition.
pm_pars <-  c("C1", "C2", "M", "P")
for (i in seq_len(nrow(pm))) {
  if (any(!(pm[i,]$est_group[[1]]$parameter %in% pm_pars))) {
    ## group_level parameters
    pm[i,]$est_group[[1]]$condition <- paste0(
      pm[i,]$est_group[[1]]$condition, 
      "_",
      str_remove(pm[i,]$est_group[[1]]$parameter, pm_pars)
    )
    rem <- str_remove(pm[i,]$est_group[[1]]$parameter, pm_pars)
    pm[i,]$est_group[[1]]$parameter <- str_remove(
      pm[i,]$est_group[[1]]$parameter, rem
    )
    
    ## individual-level parameters
    if (nrow(pm[i,]$est_indiv[[1]]) > 0) {
      pm[i,]$est_indiv[[1]]$condition <- paste0(
        pm[i,]$est_indiv[[1]]$condition, 
        "_",
        str_remove(pm[i,]$est_indiv[[1]]$parameter, pm_pars)
      )
      rem <- str_remove(pm[i,]$est_indiv[[1]]$parameter, pm_pars)
      pm[i,]$est_indiv[[1]]$parameter <- str_remove(
        pm[i,]$est_indiv[[1]]$parameter, rem
      )
    }
  }
}

check_core_pars(pm)

## add info:
pm_i <- load_combine_info(fi_pm)

all( unique(pm$dataset) %in% unique(pm_i$dataset) )
all( unique(pm$orig_model) %in% unique(pm_i$orig_model) )

pm <- left_join(pm, pm_i)
pm %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()

### convergence:
c_pm <- get_convergence(pm)
c_pm

# pm %>% 
#   group_by(dataset) %>% 
#   summarise(trials = first(trials))


##################################################################
##                          Real Model                          ##
##################################################################

## latent class missing, otherwise done
f_real <- list.files("data/real-model/", full.names = TRUE)
fi_real <- f_real[str_detect(f_real, "_info.RData$")]
f_real <- f_real[!str_detect(f_real, "_info.RData$")]
real <- load_combine(f_real)
get_info_df(real)

check_beta(real) ## all have both beta and correct prior

## remove beta cpp
real <- real %>% 
  filter(!(method == "beta++"))

## create both model variables:
real$model <- "real"
real <- real %>% 
  mutate(model2 = if_else(
    orig_model == "RealModel.eqn", "real", "real2"
  ))
real <- real %>% 
  filter(model2 == "real")

real <- check_double_fit(real)

real_core <- c("A1", "A2", "L1", "L2", "L3", "L4", "Re")
for (i in seq_len(nrow(real))) {
  ## set NA core parameters to correct value
  if (any(is.na(real[i,]$est_group[[1]]$core))) {
    real[i,]$est_group[[1]]$core <- 
      if_else(is.na(real[i,]$est_group[[1]]$core),
              if_else(
                real[i,]$est_group[[1]]$parameter %in% real_core,
                TRUE, FALSE
              ), real[i,]$est_group[[1]]$core)
  }
  if ((nrow(real[i,]$est_indiv[[1]]) > 0) &&
      any(is.na(real[i,]$est_indiv[[1]]$core))) {
    real[i,]$est_indiv[[1]]$core <- 
      if_else(is.na(real[i,]$est_indiv[[1]]$core),
              if_else(
                real[i,]$est_indiv[[1]]$parameter %in% real_core,
                TRUE, FALSE
              ), real[i,]$est_indiv[[1]]$core)
  }
}


check_core_pars(real)

## add info:
real_i <- load_combine_info(fi_real)

all( unique(real$dataset) %in% unique(real_i$dataset) )
all( unique(real$orig_model) %in% unique(real_i$orig_model) )

real <- left_join(real, real_i)
real %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()
### convergence:
c_real <- get_convergence(real)
c_real

# real %>% 
#   group_by(dataset) %>% 
#   summarise(trials = first(trials))


## check if both model variants have same core parameters
# real %>%
#   unnest(est_group) %>% 
#   group_by(model) %>% 
#   group_split() %>% 
#   map(~unique(.$parameter[.$core]))
# 
# NA_pars <- real %>%
#   unnest(est_group) %>% 
#   group_by(model, dataset, method) %>% 
#   group_split() %>% 
#   map(~unique(.$parameter[is.na(.$core)]))
# 
# real %>%
#   unnest(est_group) %>% 
#   group_by(model, dataset, method) %>% 
#   group_keys() %>% 
#   slice(which(!map_lgl(NA_pars, ~length(.) == 0)))




#################################################################
##                           r model                           ##
#################################################################

## done

f_rm <- list.files("data/r-model/", full.names = TRUE)
fi_rm <- f_rm[str_detect(f_rm, "_info.RData$")]
f_rm <- f_rm[!str_detect(f_rm, "_info.RData$")]
rm <- load_combine(f_rm)
get_info_df(rm)

check_beta(rm) ## all have correct prior

## create both model variables:
rm$model <- "rm"
rm$model2 <- "rm"

rm <- check_double_fit(rm)

c_rm <- get_convergence(rm)
c_rm

## add info:
rm_i <- load_combine_info(fi_rm)

all( unique(rm$dataset) %in% unique(rm_i$dataset) )
all( unique(rm$orig_model) %in% unique(rm_i$orig_model) )

rm <- left_join(rm, rm_i)
rm %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()

### convergence:

check_core_pars(rm)

# rm %>% 
#   group_by(dataset) %>% 
#   summarise(trials = first(trials))

##################################################################
##                        Hindsight Bias                        ##
##################################################################

## latent class impossible, because of constants

f_hb <- list.files("data/hindsight-bias/", full.names = TRUE)
f_hb <- f_hb[!str_detect(f_hb, "infos_HB.RData")]
hb <- load_combine(f_hb)
get_info_df(hb)

try(check_beta(hb)) ## fails! However, Julia says it all uses correct prior.

## remove beta++ by hand 
hb <- hb %>% 
  filter(!(method %in% c("beta++"))) 
levels(hb$method)

hb_core <- c("b", "c", "rc", "re")
for (i in seq_len(nrow(hb))) {
  ## set NA core parameters to correct value
  if (any(is.na(hb[i,]$est_group[[1]]$core))) {
    hb[i,]$est_group[[1]]$core <- 
      if_else(is.na(hb[i,]$est_group[[1]]$core),
              if_else(
                hb[i,]$est_group[[1]]$parameter %in% hb_core,
                TRUE, FALSE
              ), hb[i,]$est_group[[1]]$core)
  }
  
  ## set within-subject condition to be part of condition 
  if (any(str_detect(hb[i,]$est_group[[1]]$parameter, "w$"))) {
    nc <- str_extract(hb[i,]$est_group[[1]]$parameter, "w$")
    nc <- if_else(is.na(nc), "o", nc)
    hb[i,]$est_group[[1]]$condition <- paste0(
      hb[i,]$est_group[[1]]$condition, "_", nc
    )
    hb[i,]$est_group[[1]]$parameter <- 
      str_remove(hb[i,]$est_group[[1]]$parameter, "w$")
  }
  if ((nrow(hb[i,]$est_indiv[[1]]) > 0) && 
      any(str_detect(hb[i,]$est_indiv[[1]]$parameter, "w$"))) {
    nc <- str_extract(hb[i,]$est_indiv[[1]]$parameter, "w$")
    nc <- if_else(is.na(nc), "o", nc)
    hb[i,]$est_indiv[[1]]$condition <- paste0(
      hb[i,]$est_indiv[[1]]$condition, "_", nc
    )
    hb[i,]$est_indiv[[1]]$parameter <- 
      str_remove(hb[i,]$est_indiv[[1]]$parameter, "w$")
  }
}

hb <- check_double_fit(hb)

hb$model <- "hb"
hb$model2 <- "hb"

check_core_pars(hb)

## add info:
hb_e <- new.env()
hb_f <- load("data/hindsight-bias/infos_HB.RData", envir = hb_e)
hb_i <- map_dfr(hb_f, ~make_info_df(., env = hb_e))

all( unique(hb$dataset) %in% unique(hb_i$dataset) )
all( unique(hb$orig_model) %in% unique(hb_i$orig_model) )

hb <- left_join(hb, hb_i)
hb %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()

### following email from Julia from 06/10/2019 after discussion with Edgar,
### retain only "within" condition from Erdfelder2007
for (i in seq_len(nrow(hb))) {
  if (hb[i,]$dataset == "ErdfelderBrandtBröder2007Exp1.csv") {
    for (j in seq_along(hb)) {
      if (is.data.frame(hb[i,][[j]][[1]]) && 
          ("condition" %in% colnames(hb[i,][[j]][[1]]))) {
        #print("new")
        #print(nrow(hb[i,][[j]][[1]]))
        hb[i,][[j]][[1]] <- hb[i,][[j]][[1]] %>% 
          filter(condition == "within")
        #print(nrow(hb[i,][[j]][[1]]))
      }
    }
  }
}

#################################################################
##                       Pair Clustering                       ##
#################################################################

f_pc <- list.files("data/pair-clustering/", pattern = "RData$",
                   full.names = TRUE)
fi_pc <- f_pc[str_detect(f_pc, "_info.RData$")]
f_pc <- f_pc[!str_detect(f_pc, "_info.RData$")]
pc <- load_combine(f_pc)
pc$dataset <- str_remove(pc$dataset, pattern = ".+\\/")
pc <- pc %>% 
  mutate(dataset = case_when(
    dataset == "Broeder.csv" ~ "Broeder_sixTrials.csv",
    dataset == "Francis et al 2018_Englisch Ausschließlich.csv" ~ 
      "Francis_fourTrials_english.csv",
    dataset == "Francis et al 2018_Englisch Dominant.csv" ~ 
      "Francis_fourTrials_englishDominant.csv",
    dataset == "Francis et al 2018_Spanisch Dominant.csv" ~
      "Francis_fourTrials_spanishDominant.csv",
    dataset == "matzke.csv" ~ "Matzke_fourTrials.csv",
    dataset == "Nachmittag six Trials.csv" ~ 
      "GolzErdfelder_sixTrials_afternoon.csv",
    dataset == "Riefer.csv" ~ "Riefer_sixTrials_alcoholics.csv",
    dataset == "Riefer_Schizo.csv" ~ "Riefer_sixTrials_schizos.csv",
    dataset == "Vormittag six Trials.csv" ~ 
      "GolzErdfelder_sixTrials_forenoon.csv",
    TRUE ~ "ERROR"
  )) %>% 
  mutate(orig_model = case_when(
    orig_model %in% c("Broeder.eqn", "GolzErdfelder.eqn", "Riefer.eqn") ~
      "pc_model_sixTrials.eqn", 
    orig_model %in% c("Francis.eqn", "matzke.eqn") ~
      "pc_model_fourTrials.eqn",
    TRUE ~ "ERROR"
  ))
get_info_df(pc)

check_beta(pc)

pc <- check_double_fit(pc)

c_pc <- get_convergence(pc)
c_pc

pc$model <- "pc"
pc$model2 <- "pc"

### rename parameters
for (i in seq_len(nrow(pc))) {
  if (any(str_detect(pc[i,]$est_group[[1]]$parameter, "\\d"))) {
    
    cn <- str_extract(
      pc[i,]$est_group[[1]]$parameter, "\\d$"
    )
    cn <- if_else(is.na(cn), "1", cn)
    pc[i,]$est_group[[1]]$condition <- paste(
      pc[i,]$est_group[[1]]$condition, 
      cn, 
      sep = "_"
    )
    pc[i,]$est_group[[1]]$parameter <- str_extract(
      pc[i,]$est_group[[1]]$parameter, "[[:alpha:]]"
    )
    
  }
  ## individual-level parameters
  if ( (nrow(pc[i,]$est_indiv[[1]]) > 0) && 
       any(str_detect(pc[i,]$est_indiv[[1]]$parameter, "\\d"))) {
    
    cn <- str_extract(
      pc[i,]$est_indiv[[1]]$parameter, "\\d$"
    )
    cn <- if_else(is.na(cn), "1", cn)
    pc[i,]$est_indiv[[1]]$condition <- paste(
      pc[i,]$est_indiv[[1]]$condition, 
      cn, 
      sep = "_"
    )
    pc[i,]$est_indiv[[1]]$parameter <- str_extract(
      pc[i,]$est_indiv[[1]]$parameter, "[[:alpha:]]"
    )
  }
}



## add info:
pc_i <- load_combine_info(fi_pc)

all( unique(pc$dataset) %in% unique(pc_i$dataset) )
all( unique(pc$orig_model) %in% unique(pc_i$orig_model) )

pc <- left_join(pc, pc_i)
pc %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()

### convergence:

check_core_pars(pc)

##################################################################
##                       2-High-Threshold                       ##
##################################################################

f_c2ht <- list.files("data/2ht/", full.names = TRUE)
fi_c2ht <- f_c2ht[str_detect(f_c2ht, "_info.RData$")]
f_c2ht <- f_c2ht[!str_detect(f_c2ht, "_info.RData$")]
c2ht <- load_combine(f_c2ht)
get_info_df(c2ht)

beta_c2ht <- check_beta(c2ht)

## remove beta cpp
c2ht <- c2ht %>% 
  filter(!(method == "beta++"))

c2ht <- c2ht %>% 
  filter(method == "beta") %>% 
  anti_join(filter(beta_c2ht, is.na(prior.beta) | 
                     prior.beta != "dgamma(1,.1)T(1,)")) %>% 
  bind_rows(filter(c2ht, !(method == "beta")) ) 

c2ht <- c2ht %>% 
  mutate(model = str_remove(model, "models/cr2htm_")) %>% 
  mutate(model = str_remove(model, "\\.eqn")) %>% 
  mutate(model = str_remove(model, "point")) %>% 
  mutate(orig_model = str_remove(orig_model, "models/"))

c2ht <- check_double_fit(c2ht)

c_c2ht <- get_convergence(c2ht)
c_c2ht


c2ht_models <- c(
  "6_rrest", "8_rrest"
)

c2ht <- c2ht %>% 
  filter(model %in% c2ht_models)

c2ht$model2 <- case_when(
  c2ht$model == "6_rrest" ~ "c2ht6",
  c2ht$model == "8_rrest" ~ "c2ht8",
  TRUE ~ "ERROR"
  )
c2ht$model <- "c2ht"


check_core_pars(c2ht)

## add info:
c2ht_i <- load_combine_info(fi_c2ht)

all( sort(unique(c2ht$dataset)) %in% sort(unique(c2ht_i$dataset)) )
all( unique(c2ht$orig_model) %in% unique(c2ht_i$orig_model) )

c2ht <- left_join(c2ht, c2ht_i)
c2ht %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()


##################################################################
##                          Quad Model                          ##
##################################################################

f_qm1 <- list.files("data/quad-model/", pattern = "RData$",full.names = TRUE) 
qm1 <- load_combine(f_qm1)
qm1 <- qm1 %>% 
  mutate(dataset = case_when(
    dataset == "CalanchiniEtAl2014_PI_able" ~ "CalanchiniEtAl2014_PI_ableMCMC",
    dataset == "CalanchiniEtAl2014_PI_age" ~ "CalanchiniEtAl2014_PI_ageMCMC",
    dataset == "CalanchiniEtAl2014_PI_career" ~ 
      "CalanchiniEtAl2014_PI_careerMCMC",
    dataset == "CalanchiniEtAl2014_PI_race" ~ "CalanchiniEtAl2014_PI_raceMCMC",
    dataset == "CalanchiniEtAl2014_PI_skin" ~ "CalanchiniEtAl2014_PI_skinMCMC",
    dataset == "CalanchiniEtAl2014_PI_straight" ~ 
      "CalanchiniEtAl2014_PI_straightMCMC",
    dataset == "CalanchiniEtAl2014_UCD_High_AW" ~ 
      "CalanchiniEtAl2014_1a_AsianWhite",
    dataset == "CalanchiniEtAl2014_UCD_High_BW" ~
      "CalanchiniEtAl2014_1a_BlackWhite",
    dataset == "CalanchiniEtAl2014_UCD_Low_FI" ~ 
      "CalanchiniEtAl2014_1c_FlowerInsect",
    dataset == "CalanchiniEtAl2014_UCD_Low_ST" ~ 
      "CalanchiniEtAl2014_1c_Stereotype",
    dataset == "CalanchiniEtAl2014_UCD_Mod_BW" ~ 
      "CalanchiniEtAl2014_1b_BlackWhite",
    dataset == "CalanchiniEtAl2014_UCD_Mod_FI" ~
      "CalanchiniEtAl2014_1b_FlowerInsect",
    dataset == "WrzusEtAl2017HU" ~ "WrzusEtAl_HappyUnhappy",
    dataset == "WrzusEtAl2017NL" ~ "WrzusEtAl_NumberLetter",
    TRUE ~ dataset
  ))
get_info_df(qm1)
qm <- qm1

# qm %>%
#   filter(dataset == "LuekeGibson2015age") %>%
#   filter(inter == "Trait_u PP") %>%
#   select(convergence) %>%
#   unnest(convergence) %>%
#   filter(Rhat > 1.01)

# qm %>%
#   filter(dataset == "LuekeGibson2015age") %>%
#   filter(inter == "Trait_u PP") %>%
#   select(convergence) %>%
#   unnest(convergence) %>%
#   filter(n.eff != 0 & n.eff < 5000)


get_info_df(qm)

check_beta(qm)

qm <- check_double_fit(qm)

c_qm <- get_convergence(qm)
c_qm

qm$model <- "quad"
qm$model2 <- "quad"

## set all parameters to core parameters
for (i in seq_len(nrow(qm))) {
  qm[i,"est_group"][[1]][[1]]$core <- TRUE
}

check_core_pars(qm)

## add info:
f_qmi <- list.files("data/quad-model/get_info/", pattern = "RData$",full.names = TRUE)
qm_i <- load_combine_info(f_qmi) %>% 
  mutate(dataset = str_remove(dataset, pattern = ".csv")) %>% 
  mutate(dataset = case_when(
    dataset == "UCD_High_AW" ~ "CalanchiniEtAl2014_1a_AsianWhite",
    dataset == "UCD_High_BW" ~ "CalanchiniEtAl2014_1a_BlackWhite",
    dataset == "UCD_Mod_BW" ~ "CalanchiniEtAl2014_1b_BlackWhite",
    dataset == "UCD_Mod_FI" ~ "CalanchiniEtAl2014_1b_FlowerInsect",
    dataset == "UCD_Low_FI" ~ "CalanchiniEtAl2014_1c_FlowerInsect",
    dataset == "UCD_Low_ST" ~ "CalanchiniEtAl2014_1c_Stereotype",
    dataset == "PI_able_random" ~ "CalanchiniEtAl2014_PI_ableMCMC",
    dataset == "PI_age_random" ~ "CalanchiniEtAl2014_PI_ageMCMC",
    dataset == "PI_career_random" ~ "CalanchiniEtAl2014_PI_careerMCMC",
    dataset == "PI_race_random" ~ "CalanchiniEtAl2014_PI_raceMCMC",
    dataset == "PI_skin_random" ~ "CalanchiniEtAl2014_PI_skinMCMC",
    dataset == "PI_straight_random" ~ "CalanchiniEtAl2014_PI_straightMCMC",
    dataset == "WrzusEtAl2017HU" ~ "WrzusEtAl_HappyUnhappy",
    dataset == "WrzusEtAl2017NL" ~ "WrzusEtAl_NumberLetter",
    TRUE ~ dataset
  ))

all( sort(unique(qm$dataset)) %in% sort(unique(qm_i$dataset)) )
all( sort(unique(qm_i$dataset)) %in% sort(unique(qm$dataset)) )
all( unique(qm$orig_model) %in% unique(qm_i$orig_model) )

qm <- left_join(qm, qm_i)

#################################################################
##                    2HTSM - Source Memory                    ##
#################################################################


## submodel 4
f_sm4 <- list.files("data/2htsm/Submodel 4/", 
                    pattern = "RData$",full.names = TRUE)
htsm_4 <- load_combine(f_sm4)
htsm_4$dataset <- str_remove_all(htsm_4$dataset, "\\.csv|\\.cvs")
htsm_4 <- htsm_4 %>% 
  mutate(dataset = case_when(
    dataset == "K2016" ~ "K2016_E1",
    dataset == "KT017" ~ "KT2017",
    TRUE ~ dataset
  )) 

get_info_df(htsm_4)

check_beta(htsm_4) 

htsm_4 <- htsm_4 %>% 
  filter(!(method == "beta++"))


htsm_4 %>%
  filter(dataset == "MH2001") %>%
  unnest(est_group) %>%
  group_by(inter, orig_condition) %>%
  tally %>% 
  print(n = Inf)
  
htsm_4 <- check_double_fit(htsm_4)
htsm_4$model <- "2htsm"
htsm_4$model2 <- "2htsm_4"

get_convergence(htsm_4)
check_core_pars(htsm_4)

## submodel 5d
f_sm5d <- list.files("data/2htsm/Submodel 5d/", 
                    pattern = "RData$",full.names = TRUE)
htsm_5d <- load_combine(f_sm5d)
htsm_5d$dataset <- str_remove_all(htsm_5d$dataset, "\\.csv|\\.cvs")
get_info_df(htsm_5d)

check_beta(htsm_5d) %>% 
  filter(is.na(prior.beta) | prior.beta != "dgamma(1,.1)T(1,)") %>% 
  ungroup() %>% 
  summarise(n_beta_wrong_prior = sum(n_beta))

htsm_5d <- htsm_5d %>% 
  filter(!(method == "beta++"))

htsm_5d <- check_double_fit(htsm_5d)
htsm_5d$model <- "2htsm"
htsm_5d$model2 <- "2htsm_5d"

get_convergence(htsm_5d)

check_core_pars(htsm_5d)

## submodel 6e
f_sm6e <- list.files("data/2htsm/Submodel 6e/", 
                    pattern = "RData$",full.names = TRUE)
htsm_6e <- load_combine(f_sm6e) %>% 
  mutate(model = "2htsm", 
         model2 = "2htsm_6e") %>% 
  mutate(orig_model = case_when(
    orig_model == "2HTSM_Sub6e_2WSCond.eqn" ~ "2HTSM_Buchnermodel_2WSCond.eqn",
    TRUE ~ orig_model
  ))
htsm_6e$dataset <- str_remove_all(htsm_6e$dataset, "\\.csv|\\.cvs")
get_info_df(htsm_6e)

htsm_6e$orig_model %>% 
  unique

check_beta(htsm_6e) %>%
  filter(is.na(prior.beta) | prior.beta != "dgamma(1,.1)T(1,)") %>% 
  ungroup() %>% 
  summarise(n_beta_wrong_prior = sum(n_beta))

htsm_6e <- htsm_6e %>% 
  filter(!(method == "beta++"))

htsm_6e <- check_double_fit(htsm_6e)

get_convergence(htsm_6e)

check_core_pars(htsm_6e)

### combine

htsm <- bind_rows(
  htsm_4, htsm_5d, htsm_6e
)

### latent class (removed after general agreement on the matter)

# f_htsm_lc <- list.files("data/2htsm/latent class/", 
#                     pattern = "RData$",full.names = TRUE)
# htsm_lc <- load_combine(f_htsm_lc)
# htsm_lc$dataset <- str_remove_all(htsm_lc$dataset, "\\.csv|\\.cvs")
# 
# ## remove model 5a and rename 6e to Buchnermodel
# htsm_lc <- htsm_lc %>% 
#   filter(!(model %in% c("2HTSM_Sub5a_1WSCond.eqn"))) %>% 
#   mutate(orig_model = case_when(
#     orig_model == "2HTSM_Sub6e_1WSCond.eqn" ~ "2HTSM_Buchnermodel_1WSCond.eqn",
#     TRUE ~ orig_model
#   ))
# 
# 
# # htsm_lc %>% 
# #   filter(dataset == "BK2011_E1") %>% 
# #   unnest(est_group) %>% 
# #   group_by(inter, orig_condition) %>% 
# #   tally
# 
# htsm_lc %>% 
#   group_by(dataset, inter) %>% 
#   tally %>% 
#   spread(inter, n) %>% 
#   print(n = Inf)
# 
# 
# all(unique(htsm_lc$dataset) %in% unique(htsm$dataset)) ## should be TRUE
# 
# htsm_lc <- htsm_lc %>% 
#   mutate(model = "2htsm", 
#          model_type = case_when(
#            dataset %in% unique(htsm_4$dataset) ~ "main",
#            TRUE ~ "sub"
#          ))


### combine

## exclude latent class due to inconsistent naming of condition
htsm <- bind_rows(
  htsm_4, htsm_5d, htsm_6e #, htsm_lc
)

htsm <- check_double_fit(htsm)

## set all parameters to core parameters and move within-subject 
## note: last number always corresponds to condition!
### number to condition
for (i in seq_len(nrow(htsm))) {
  htsm[i,"est_group"][[1]][[1]]$core <- TRUE
    if (htsm[i,]$dataset == "BK2011_E2") {
    for (j in seq_along(htsm)) {
      if (is.data.frame(htsm[i,j][[1]][[1]])) {
        if ("orig_condition" %in% colnames(htsm[i,j][[1]][[1]])) {
          htsm[i,j][[1]][[1]]$orig_condition <- case_when(
            htsm[i,j][[1]][[1]]$orig_condition == "1" ~ "cons",
            htsm[i,j][[1]][[1]]$orig_condition == "2" ~ "incons",
            htsm[i,j][[1]][[1]]$orig_condition == "3" ~ "zero",
            TRUE ~ htsm[i,j][[1]][[1]]$orig_condition
          )
        }
        if ("condition" %in% colnames(htsm[i,j][[1]][[1]])) {
          htsm[i,j][[1]][[1]]$condition <- case_when(
            htsm[i,j][[1]][[1]]$condition == "1" ~ "cons",
            htsm[i,j][[1]][[1]]$condition == "2" ~ "incons",
            htsm[i,j][[1]][[1]]$condition == "3" ~ "zero",
            TRUE ~ htsm[i,j][[1]][[1]]$condition
          )
        }
      } 
    }
  }
  if (htsm[i,]$dataset == "MH2001") {
    for (j in seq_along(htsm)) {
      if (is.data.frame(htsm[i,j][[1]][[1]])) {
        if ("orig_condition" %in% colnames(htsm[i,j][[1]][[1]])) {
          htsm[i,j][[1]][[1]]$orig_condition <- case_when(
            htsm[i,j][[1]][[1]]$orig_condition == "2" ~ "group_sources",
            htsm[i,j][[1]][[1]]$orig_condition == "3" ~ "town_sources",
            TRUE ~ htsm[i,j][[1]][[1]]$orig_condition
          )
        }
        if ("condition" %in% colnames(htsm[i,j][[1]][[1]])) {
          htsm[i,j][[1]][[1]]$condition <- case_when(
            htsm[i,j][[1]][[1]]$condition == "2" ~ "group_sources",
            htsm[i,j][[1]][[1]]$condition == "3" ~ "town_sources",
            TRUE ~ htsm[i,j][[1]][[1]]$condition
          )
        }
      } 
    }
  }
  if (any(str_detect(htsm[i,]$est_group[[1]]$parameter, "\\d"))) {
    if (str_detect(htsm[i,]$orig_model, "_Sub4_")) {
      cn <- str_extract(htsm[i,]$est_group[[1]]$parameter, "\\d")
      cn <- if_else(is.na(cn), "1", cn)
      htsm[i,]$est_group[[1]]$condition <- paste(
        htsm[i,]$est_group[[1]]$condition, 
        cn, 
        sep = "_"
      )
      htsm[i,]$est_group[[1]]$parameter <- 
        str_remove(htsm[i,]$est_group[[1]]$parameter, "\\d")
    } else {
      cn <- str_extract(htsm[i,]$est_group[[1]]$parameter, "[[:alpha:]]\\d") %>% 
        substr(2, 2)
      cn <- if_else(is.na(cn), "1", cn)
      htsm[i,]$est_group[[1]]$condition <- paste(
        htsm[i,]$est_group[[1]]$condition, 
        cn, 
        sep = "_"
      )
      htsm[i,]$est_group[[1]]$parameter <- if_else(
        str_detect(htsm[i,]$est_group[[1]]$parameter, "[[:alpha:]]\\d"),
        str_remove(htsm[i,]$est_group[[1]]$parameter, "\\d"),
        htsm[i,]$est_group[[1]]$parameter
        
      ) 
      
    }
  }
  
  ## individual-level parameters
  if ( (nrow(htsm[i,]$est_indiv[[1]]) > 0) && 
       any(str_detect(htsm[i,]$est_indiv[[1]]$parameter, "\\d"))) {
    if (str_detect(htsm[i,]$orig_model, "_Sub4_")) {
      cn <- str_extract(htsm[i,]$est_indiv[[1]]$parameter, "\\d")
      cn <- if_else(is.na(cn), "1", cn)
      htsm[i,]$est_indiv[[1]]$condition <- paste(
        htsm[i,]$est_indiv[[1]]$condition, 
        cn, 
        sep = "_"
      )
      htsm[i,]$est_indiv[[1]]$parameter <- 
        str_remove(htsm[i,]$est_indiv[[1]]$parameter, "\\d")
    } else {
      cn <- str_extract(htsm[i,]$est_indiv[[1]]$parameter, "[[:alpha:]]\\d") %>% 
        substr(2, 2)
      cn <- if_else(is.na(cn), "1", cn)
      htsm[i,]$est_indiv[[1]]$condition <- paste(
        htsm[i,]$est_indiv[[1]]$condition, 
        cn, 
        sep = "_"
      )
      htsm[i,]$est_indiv[[1]]$parameter <- if_else(
        str_detect(htsm[i,]$est_indiv[[1]]$parameter, "[[:alpha:]]\\d"),
        str_remove(htsm[i,]$est_indiv[[1]]$parameter, "\\d"),
        htsm[i,]$est_indiv[[1]]$parameter
      ) 
      
    }
  }
}

## change d_3 into d_2 in DS2000_E1:
for (i in which(htsm$dataset == "DS2000_E1")) {
  ## est group:
  htsm[i,"est_group"][[1]][[1]]$parameter <- if_else(
    htsm[i,"est_group"][[1]][[1]]$parameter == "d_3", 
    "d_2", 
    htsm[i,"est_group"][[1]][[1]]$parameter
  )
  # also set orig_parameter for consistency to d3_3
  htsm[i,"est_group"][[1]][[1]]$orig_parameter <- if_else(
    htsm[i,"est_group"][[1]][[1]]$orig_parameter == "d3_3",
    "d3_2",
    htsm[i,"est_group"][[1]][[1]]$orig_parameter
  )
  ## est indiv:
  if ( (nrow(htsm[i,]$est_indiv[[1]]) > 0)) {
    htsm[i,"est_indiv"][[1]][[1]]$parameter <- if_else(
      htsm[i,"est_indiv"][[1]][[1]]$parameter == "d_3", 
      "d_2", 
      htsm[i,"est_indiv"][[1]][[1]]$parameter
    )
    # also set orig_parameter for consistency to d3_3
    htsm[i,"est_indiv"][[1]][[1]]$orig_parameter <- if_else(
      htsm[i,"est_indiv"][[1]][[1]]$orig_parameter == "d3_3",
      "d3_2",
      htsm[i,"est_indiv"][[1]][[1]]$orig_parameter
    )
  }
}

check_core_pars(htsm)

# htsm %>%
#   unnest(est_group) %>%
#   select(parameter, orig_parameter) %>%
#   unique %>%
#   print(n = Inf)
# 
# htsm %>%
#   unnest(est_group) %>% 
#   filter(orig_parameter == "d3_3") %>% 
#   select(model, dataset, condition, parameter, orig_parameter)

## add info:
htsm_e <- new.env()
htsm_f <- load("data/2htsm/infos_2HTSM.RData", envir = htsm_e)
htsm_f <- htsm_f[-8]
htsm_i <- map_dfr(htsm_f, ~make_info_df(., env = htsm_e))

htsm_i <- htsm_i %>% 
  mutate(dataset = str_remove(dataset, ".csv$")) %>% 
  mutate(dataset = str_replace_all(dataset, "ö", "oe")) %>% 
  mutate(dataset = str_replace_all(dataset, "ü", "ue")) %>% 
  mutate(dataset = if_else(dataset == "Schuetz & Broeder (2011) Ex2", 
                           "Schuetz & Broeder", dataset)) %>% 
  mutate(dataset = if_else(dataset == "Kueppers & Bayen (2014)", 
                           "Kueppers & Bayen 2014", dataset))

all( sort(unique(htsm$dataset)) %in% sort(unique(htsm_i$dataset)) )

all( sort(unique(htsm$orig_model)) %in% sort(unique(htsm_i$orig_model)) )

htsm <- left_join(htsm, htsm_i)
htsm %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()

htsm %>% 
  filter(orig_model == "2HTSM_Sub6e_1WSCond.eqn") %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree)

c_htsm <- get_convergence(htsm)
c_htsm


## check to see if each data set is only fitted using one model type:
htsm %>% 
  select(dataset, model2) %>% 
  unique %>% 
  table() %>% 
  apply(1, function(x) sum(x) > 1) %>% 
  any()
## FALSE indicates each data set is fitted using one model type


htsm %>% 
  filter(dataset == "BK2011_E2") %>% 
  unnest(est_group) %>%
  group_by(inter, orig_condition) %>%
  tally

htsm %>% 
  unnest(est_group) %>%
  summarise(any(is.na(condition)))

# htsm %>% 
#   filter(dataset == "DS2000_E1") %>% 
#   select(model_exp) %>% 
#   slice(1) %>% 
#   {.[[1]]}

##################################################################
##                     Process Dissociation                     ##
##################################################################

f_pd <- list.files("data/process-dissociation/", full.names = TRUE)
fi_pd <- f_pd[str_detect(f_pd, "_info.RData$")]
f_pd <- f_pd[!str_detect(f_pd, "_info.RData$")]

load("data/process-dissociation/for-henrik.RData")
pd <- sanitize_vars(prepare_single_result(results))
pd <- pd %>% 
  filter(!(method == "beta++"))

pd$orig_model <- pd$model
rename_models <- c(
  bodner_2000_exp2.eqn = "pd_s"
  , bodner_2000_exp4.eqn = "pd_s"
  , caldwell_masson_2001_exp1.eqn = "pd_s"
  , caldwell_masson_2001_exp2.eqn = "pd_s"
  , `pd-alternative.eqn` = "pd_e"
  , klauer_2015_exp3.eqn = "pd_e"
  , klauer_2015_exp5_m1.eqn = "pd_e"
  , rouder_2008_exp1.eqn = "pd_s"
  , rouder_2008_exp2.eqn = "pd_s"
  , stahl_2015_exp1.eqn = "pd_s"
)
pd$model2 <- rename_models[pd$model]
pd$model <- "pd"
# table(pd$model2)

new_par_names <- c(
  a1 = "A" #1
  , a2 = "A" #2
  , a3 = "A" #3
  , c1 = "C" #1
  , c2 = "C" #2
  , c3 = "C" #3
  , c = "C" #1
  , u = "A" #1
  , a = "A" #1
  , H_50 = "A"
  , H_75 = "A"
  , R_50 = "C"
  , R_75 = "C"
  , h1 = "A_alt"
  , h2 = "A_alt"
  , h3 = "A_alt"
  , wi = "C_Exclusion"
  , wc = "C_Inclusion"
  , A = "A"
  , C = "C"
  , ab1 = "A_alt" # ws1
  , ab2 = "A_alt" # ws2
  , ab3 = "A_alt" # ws3
  , aw1 = "A_alt" # ws4
  , aw2 = "A_alt" # ws5
  , aw3 = "A_alt" # ws6
  , bg = "C_Inclusion"
  , bt = "C_Exclusion"
  , wg = "C_Exclusion"
  , wt = "C_Inclusion"
)

new_cond_names <- c(
    # Bodner et al. ----
  a1 = 1
  , a2 = 2
  , a3 = 3
  , c1 = 1
  , c2 = 2
  , c3 = 3
  , c = 1
  , u = 1
  , a = 1
  , H_50 = 1
  , H_75 = 2
  , R_50 = 1
  , R_75 = 2
  , A = 1
  , C = 1
  # Klauer et al. ----
  , h1 = 1
  , h2 = 2
  , h3 = 3
  , wi = 1
  , wc = 1
  , ab1 = 1
  , ab2 = 2
  , ab3 = 3
  , aw1 = 4
  , aw2 = 5
  , aw3 = 6
  , bg = 1
  , bt = 1
  , wg = 2
  , wt = 2
)


for (i in seq_len(nrow(pd))) { 
  pd[i,]$est_group[[1]]$parameter <- new_par_names[pd[i,]$est_group[[1]]$orig_parameter]
  pd[i,]$est_group[[1]]$condition <- paste0(
    pd[i,]$est_group[[1]]$orig_condition, "_", 
    new_cond_names[pd[i,]$est_group[[1]]$orig_parameter]
  )
  if ( (nrow(pd[i,]$est_indiv[[1]]) > 0) ) {
    pd[i,]$est_indiv[[1]]$parameter <- new_par_names[pd[i,]$est_indiv[[1]]$parameter]
    pd[i,]$est_indiv[[1]]$condition <- paste0(
      pd[i,]$est_indiv[[1]]$orig_condition, "_", 
      new_cond_names[pd[i,]$est_indiv[[1]]$orig_parameter]
    )
  }
}

rm(results)
# pd <- load_combine(f_pd)

get_info_df(pd)
check_beta(pd)

pd <- check_double_fit(pd)

## add info:
pd_i <- load_combine_info(fi_pd)

all( sort(unique(pd$dataset)) %in% sort(unique(pd_i$dataset)) )
all( unique(pd$orig_model) %in% unique(pd_i$orig_model) )

pd <- left_join(pd, pd_i)
pd %>% 
  select(n_tree, rel_tree, model_df, model_exp, data_tree) %>% 
  summarise_all(list(null = ~any(map_lgl(., is.null)), 
                     na = ~any(map_lgl(., ~isTRUE(is.na(.)))))) %>% 
  as.data.frame()

check_core_pars(pd)


##----------------------------------------------------------------
##                          Combine Data                         -
##----------------------------------------------------------------

dm <- bind_rows(
  pm,
  real, 
  rm,
  hb,
  htsm, 
  c2ht,
  pc
  , pd
  , qm
)


### trials are missing for most data sets

# c_dm <- left_join(c_pm, c_real) %>% 
#   left_join(c_rm) %>% 
#   left_join(c_htsm) %>% 
#   left_join(c_c2ht)

new_levels <- c("Comp MLE", "Comp Bayes", 
                                 "No asy", "No PB", "No NPB", "No Bayes",  
                                 "LC PP",
                                 "Beta PP", 
                                 "Trait_u PP", "Trait PP")

dm <- dm %>% 
  mutate(inter = factor(inter, levels = new_levels))

dm %>% 
  filter(is.na(inter)) %>% 
  print(n = 100)

dput(unique(dm$model2))
any(is.na(dm$model2))

dput(unique(dm$model))

model_levels <- c("2htsm", "c2ht", "pc", "pd", "pm", "hb", "rm", "real", "quad")
model2_levels <- c("2htsm_4", "2htsm_5d", "2htsm_6e", 
                   "c2ht6", "c2ht8", "pc", "pd_s", "pd_e", 
                   "pm", "hb", "rm", "real", "quad")

dm <- dm %>% 
  mutate(model2 = factor(model2, levels = model2_levels),
         model = factor(model, levels = model_levels)
  )
stopifnot(!any(is.na(dm$model)))  ## must be TRUE
stopifnot(!any(is.na(dm$model2)))

dm %>% 
  group_by(model) %>% 
  tally

# dm %>% 
#   filter(is.na(model2)) %>% 
#   group_by(model) %>% 
#   tally

# dm_l <- list(
#   pm = pm,
#   rm = rm,
#   hb = hb,
#   htsm = htsm, 
#   c2ht = c2ht,
#   pc = pc,
#   real = real
#   #, quad = qm
# )
# 
# dm_l <- map(dm_l, ~{
#   .$inter <- factor(.$inter, levels = new_levels)
#   .$model_type = factor(.$model_type, levels = c("main", "sub"))
#   .$model = factor(.$model, levels = model_levels)
#   .
# })

###

n_participants <- dm %>%
  unnest(est_indiv) %>%
  group_by(model, dataset, orig_condition, inter, model2, id) %>%
  summarise(parameter = first(parameter)) %>%
  tally(name = "n_participant")

n_participants %>%
  group_by(model, dataset, orig_condition, model2) %>%
  summarise(check = all(n_participant == first(n_participant))) %>%
  {all(.$check)}  %>% 
  stopifnot

# n_participants %>%
#   group_by(model, dataset, orig_condition, model2) %>%
#   summarise(check = all(n_participant == first(n_participant))) %>% 
#   filter(!check)

# n_participants %>%
#   group_by(model, dataset, orig_condition, model2) %>%
#   filter(dataset == "CalanchiniEtAl2014_PI_skinMCMC")
# 
# n_participants %>%
#   group_by(model, dataset, orig_condition, model2) %>%
#   filter(dataset == "CalanchiniEtAl2014_PI_straightMCMC")

n_participants <- n_participants %>%
  group_by(model, dataset, orig_condition, model2) %>%
  summarise(n_participant = first(n_participant)) %>%
  ungroup()
# 
# n_study <- dm %>% 
#   group_by(model) %>% 
#   summarise(n_study = length(unique(dataset)))

# dm %>% 
#   unnest(gof_indiv) %>% 
#   print(n = 100)

save(dm, n_participants, 
     file = "combined_results.RData", compress = "xz")  #

