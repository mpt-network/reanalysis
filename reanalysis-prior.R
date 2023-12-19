library("TreeBUGS")
library("tidyverse")

data("roc6", package = "MPTinR")

unique(roc6$exp)

dj <- roc6 %>% 
  filter(exp == "Jaeger_2012")

get_means_2ht <- function(model) {
  sum <- summary(model)
  sum$groupParameters$mean[1:3,"Mean"]
}

res_orig <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                      mu = "dnorm(0,1)", xi = "dunif(0,10)", 
                     n.iter = 250000, n.adapt = 50000, n.burnin = 30000, n.thin = 20)

res_orig_low <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                     mu = "dnorm(-1,1)", xi = "dunif(0,10)", 
                     n.iter = 250000, n.adapt = 50000, n.burnin = 30000, n.thin = 20)

res_orig_high <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                         mu = "dnorm(1,1)", xi = "dunif(0,10)", 
                         n.iter = 250000, n.adapt = 50000, n.burnin = 30000, n.thin = 20)

cbind(
  orig = get_means_2ht(res_orig),
  orig_low = get_means_2ht(res_orig_low),
  orig_high = get_means_2ht(res_orig_high)
)

res_narrow <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                      mu = "dnorm(0,0.33)", xi = "dunif(0,10)", 
                      n.iter = 150000, n.adapt = 40000, n.burnin = 30000, n.thin = 10)
summary(res_narrow) ## ok

res_narrow_low <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                       mu = "dnorm(-1,0.33)", xi = "dunif(0,10)", 
                       n.iter = 500000, n.adapt = 40000, n.burnin = 30000, n.thin = 10)
summary(res_narrow_low) ## ok

res_narrow_high <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                           mu = "dnorm(1,0.33)", xi = "dunif(0,10)", 
                           n.iter = 500000, n.adapt = 40000, n.burnin = 30000, n.thin = 10)
summary(res_narrow_high) ## ok

res_wide <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                      mu = "dnorm(0,3)", xi = "dunif(0,10)", 
                     n.iter = 250000, n.adapt = 40000, n.burnin = 30000, n.thin = 20)
summary(res_wide)  ## ok

res_wide_low <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                     mu = "dnorm(-1,3)", xi = "dunif(0,10)", 
                     n.iter = 250000, n.adapt = 40000, n.burnin = 30000, n.thin = 20)
summary(res_wide_low) ## ok

res_wide_high <- traitMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                     mu = "dnorm(1,3)", xi = "dunif(0,10)", 
                     n.iter = 500000, n.adapt = 40000, n.burnin = 30000, n.thin = 10)
summary(res_wide_high) ## ok

prio_pars_jaeger <- cbind(
  orig = get_means_2ht(res_orig),
  orig_low = get_means_2ht(res_orig_low),
  orig_high = get_means_2ht(res_orig_high),
  narrow = get_means_2ht(res_narrow),
  narrow_low = get_means_2ht(res_narrow_low),
  narrow_high = get_means_2ht(res_narrow_high),
  wide = get_means_2ht(res_wide),
  wide_low = get_means_2ht(res_wide_low),
  wide_high = get_means_2ht(res_wide_high)
)

mean(abs(prio_pars_jaeger[,-1] - prio_pars_jaeger[,1]))
max(abs(prio_pars_jaeger[,-1] - prio_pars_jaeger[,1]))

### beta
res_beta <- betaMPT("model-discussion/cr2htm_6point_rrest.eqn", data = dj, 
                     n.iter = 250000, n.adapt = 50000, n.burnin = 30000, n.thin = 10)
summary(res_beta) ## ok

##### Baayen & Kuhlmann (2011), Exp. 1

db <- read_csv("model-discussion/BK2011_E1.csv")

resb_orig <- res_orig <- traitMPT("model-discussion/2HTSM_Sub4_1WSCond.eqn", data = db, 
                                  mu = "dnorm(0,1)", xi = "dunif(0,10)", 
                                  n.iter = 250000, n.adapt = 50000, n.burnin = 30000, n.thin = 20)
summary(resb_orig) ## ok
bk_orig <- summary(resb_orig)$groupParameters$mean[,"Mean"]

priors <- c("dnorm(-1,1)", "dnorm(1,1)", 
            "dnorm(0,0.33)", "dnorm(-1,0.33)", "dnorm(1,0.33)", 
            "dnorm(0,3)", "dnorm(-1,3)", "dnorm(1,3)")
bk_out <- vector("list", length(priors))
for (i in seq_along(bk_out)) {
  bk_out[[i]] <-  traitMPT("model-discussion/2HTSM_Sub4_1WSCond.eqn", data = db, 
                           mu = priors[i], xi = "dunif(0,10)", 
                           n.iter = 250000, n.adapt = 75000, n.burnin = 30000, n.thin = 25)
}
any(map_lgl(bk_out, ~any(summary(.)$groupParameters$mean[,"Rhat"] > 1.05))) ## FALSE means all have converged

bk_means <- map(bk_out, ~summary(.)$groupParameters$mean[,"Mean"]) 
names(bk_means) <- priors
bk_means <- bind_cols(bk_means)
bk_means - summary(resb_orig)

##### Coolin et al. (2016)

dc <- read_csv2("model-discussion/Coolin2016.csv")

resc_orig <- traitMPT("model-discussion/HB_Coolin2016.eqn", data = dc, 
                                  mu = "dnorm(0,1)", xi = "dunif(0,10)", 
                                  n.iter = 500000, n.adapt = 50000, n.burnin = 30000, n.thin = 25)
summary(resc_orig) ## 
cool_orig <- summary(resc_orig)$groupParameters$mean[,"Mean"]

priors <- c("dnorm(-1,1)", "dnorm(1,1)", 
            "dnorm(0,0.33)", "dnorm(-1,0.33)", "dnorm(1,0.33)", 
            "dnorm(0,3)", "dnorm(-1,3)", "dnorm(1,3)")
cool_out <- vector("list", length(priors))
for (i in seq_along(cool_out)) {
  cool_out[[i]] <-  traitMPT("model-discussion/HB_Coolin2016.eqn", data = dc, 
                           mu = priors[i], xi = "dunif(0,10)", 
                           n.iter = 500000, n.adapt = 75000, n.burnin = 30000, n.thin = 25)
}
any(map_lgl(cool_out, ~any(summary(.)$groupParameters$mean[,"Rhat"] > 1.05))) ## FALSE means all have converged

cool_means <- map(cool_out, ~summary(.)$groupParameters$mean[,"Mean"]) 
names(cool_means) <- priors
cool_means <- bind_cols(cool_means)
cool_means - summary(resc_orig)

save(cool_means, cool_orig, 
     bk_orig, bk_means, 
     prio_pars_jaeger, file = "prior_recalc.rda")
