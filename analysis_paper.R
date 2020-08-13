
library("afex")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))

load("all_pairs_core.RData")
all_pairs <- all_pairs %>% 
  mutate(rel_weight = (rel_par_weight_x + rel_par_weight_y)/2,
         rel_n = (rel_n_x + rel_n_y) / 2) %>% 
  mutate(se_x_w = if_else(se_x > 0.25, 0.25, se_x),
         se_y_w = if_else(se_y > 0.25, 0.25, se_y))


##----------------------------------------------------------------
##              Correlation of Relevant Covariates               -
##----------------------------------------------------------------


str(all_pairs)
targ_cmle <- all_pairs %>% 
  filter(cond_y == "Comp MLE") %>% 
  filter(cond_x != "Comp MLE")

rel_vars <- c("se_x_w", "se_y_w", "fungi_max", "rel_weight", "rel_n")

targ_cmle %>% 
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
  
