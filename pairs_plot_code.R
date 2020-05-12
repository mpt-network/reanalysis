
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
                      levels(all_pars$inter)[unique(all_pars2$cond_x) %in% 
                                               levels(all_pars$inter)]),
    cond_y = factor(cond_y, levels = 
                      levels(all_pars$inter)[unique(all_pars2$cond_y) %in% 
                                               levels(all_pars$inter)])
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
