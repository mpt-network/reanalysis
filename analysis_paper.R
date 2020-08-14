
library("afex")
library("DescTools")
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))

load("all_pairs_core.RData")
all_pairs <- all_pairs %>% 
  mutate(rel_weight = (rel_par_weight_x + rel_par_weight_y)/2,
         rel_n = (rel_n_x + rel_n_y) / 2) %>% 
  mutate(se_x_w = if_else(se_x > 0.25, 0.25, se_x),
         se_y_w = if_else(se_y > 0.25, 0.25, se_y)) %>% 
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
##              Correlation of Relevant Covariates               -
##----------------------------------------------------------------


str(all_pairs)
targ_cmle <- all_pairs %>% 
  filter(cond_y == "Comp MLE") %>% 
  filter(cond_x != "Comp MLE")

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
  
