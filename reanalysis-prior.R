library("TreeBUGS")
library("tidyverse")

data("roc6", package = "MPTinR")

unique(roc6$exp)

dj <- roc6 %>% 
  filter(exp == "Jaeger_2012")

res_orig <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                      mu = "dnorm(0,1)", xi = "dunif(0,10)", 
                     n.iter = 75000, n.adapt = 20000, n.burnin = 30000, n.thin = 10)
summary(res_orig)

res_narrow <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                      mu = "dnorm(0,0.2)", xi = "dunif(0,5)", 
                      n.iter = 100000, n.adapt = 40000, n.burnin = 30000, n.thin = 10)

res_wide <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                      mu = "dnorm(0,2.5)", xi = "dunif(0,10)", 
                     n.iter = 100000, n.adapt = 40000, n.burnin = 30000, n.thin = 10)
