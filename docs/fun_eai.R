make_emp_biv_plot <- function(data, var1, label1, var2, label2) {
  sum1 <- data %>% 
    filter(is.finite(!!enquo(var1)), is.finite(!!enquo(var2))) %>% 
    mutate(across(c({{var1}}), 
                  .fns = ~cut_width(., width = CUT_WIDTH_RHO, boundary = 0))) %>% 
    mutate(across(c({{var2}}), 
                  .fns = ~cut_width(., width = CUT_WIDTH_SE, boundary = 0))) %>% 
    group_by(cond_label, cond_iv_label, {{var2}}, {{var1}}) %>% 
    summarise(mean_abs_dev = mean(abs_dev),
              max_abs_dev = max(abs_dev),
              n = n(), 
              sd = sd(abs_dev),
              se = sd(abs_dev)/sqrt(n())) 
  
  var1name <- all.vars(enquo(var1))
  var2name <- all.vars(enquo(var2))
  sum1[[var1name]] <- make_cut_num(sum1[[var1name]])
  sum1[[var2name]] <- make_cut_num(sum1[[var2name]])
  
  sum1 %>% 
    mutate(mean_abs_dev = if_else(n >= 3, mean_abs_dev, NA_real_)) %>% 
    ggplot(aes(x = {{var2}}, y = {{var1}}, fill = mean_abs_dev)) +
    geom_raster() +
    scale_fill_gradientn(limits = c(0, 0.3), 
                         colours = c("darkgreen", "yellow", 
                                     "orange", "red", "darkred"), 
                         values = 
                           scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2, 0.3)),
                         name = "Mean abs. deviation", 
                         na.value = "transparent") +
    facet_grid(~cond_iv_label) +
    coord_fixed(ratio = 1, xlim = c(0, 0.5), ylim = c(0, 0.9)) +
    scale_x_continuous(guide = guide_axis(angle = -90, check.overlap = TRUE)) +
    # scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
    labs(x = label2, y = label1) 
}

make_gam_biv_plot <- function(data, var1, label1, var2, label2) {
  
  bigams <- data %>% 
    filter(abs_dev != 0) %>% 
    group_by(cond_label, cond_iv_label) %>% 
    nest() %>% 
    mutate(
      res = map(data, 
                ~mgcv::gam(formula(expr(
                  abs_dev ~ te(!!enexpr(var2), !!enexpr(var1), bs = "ts")
                  )), data = ., family=Gamma(link=log))
      ))
  
  var1name <- all.vars(enquo(var1))
  var2name <- all.vars(enquo(var2))
  te_term <- paste0("te(", var2name, ",",var1name, ")")
  
  bigams <- bigams %>% 
    mutate(nd = map(res, ~gratia::smooth_estimates(., te_term, n = 100))) %>% 
    mutate(nd = map2(res, nd, ~mutate(.y, abs_dev = predict(.x, newdata = .y, 
                                                            type = "response"))))
  ## mgcv::gam(formula(expr(abs_dev ~ s(!!enexpr(covariate), bs = "ts"))
  pout <- bigams %>% 
    select(cond_label, cond_iv_label, nd) %>% 
    unnest(nd) %>% 
    mutate(abs_dev = if_else(is.na(est), NA_real_, abs_dev)) %>% 
    ggplot(., aes(x = {{var2}}, y = {{var1}})) +
    geom_raster(mapping = aes(fill = abs_dev)) +
    geom_contour(aes(z = abs_dev), 
                 breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3)) +
    scale_fill_gradientn(limits = c(0, 0.3), 
                         colours = c("darkgreen", "yellow", 
                                     "orange", "red", "darkred"), 
                         values = 
                           scales::rescale(c(0, 0.01, 0.05, 0.1, 0.2, 0.3)),
                         name = "Mean abs. deviation", 
                         na.value = "transparent") +
    coord_fixed(ratio = 1, xlim = c(0, 0.5), ylim = c(0, 0.9)) +
    labs(x = label2, y = label1)  +
    scale_x_continuous(guide = guide_axis(angle = -90, check.overlap = TRUE))
  pout
}

make_cor_plot <- function(data, x, y, filter, alpha = 0.3) {
  orig_n <- nrow(data)
  if (!missing(filter)) {
    data <- data %>% 
      dplyr::filter( !!enquo(filter) )
  }
  
  varname_x <- all.vars(enquo(x))
  varname_y <- all.vars(enquo(y))
  
  data <- data %>% 
    dplyr::filter(!is.na(!!sym(varname_x)), !is.na(!!sym(varname_y)))
  
  removed <- 1 - nrow(data)/orig_n
  
  dlm <- data %>% 
    group_by(cond_iv_label) %>% 
    summarise(lm = list(lm(formula(expr(!!enexpr(y) ~ !!enexpr(x))))),
              r = cor({{y}}, {{x}})) 
  dlm <- dlm %>% 
    mutate(est = map(lm, ~ broom::tidy(.)), 
           gof = map(lm, ~ broom::glance(.)))
  dgof <- dlm %>% 
    select(-lm, -est) %>% 
    unnest(gof)
  rx <- range(data[[varname_x]], na.rm = TRUE)
  ry <- range(data[[varname_y]], na.rm = TRUE)
  
  new_xlab <- if (varname_x == as_label(enquo(x))) {
    varname_x
  } else {
    paste0(varname_x, " (", as_label(enquo(x)), ")")
  }
  if (!missing(filter)) {
    new_xlab <- paste0(new_xlab, "; filter: ", as_label(enquo(filter)))
  }
  if (removed > 0) {
    new_xlab <- paste0(new_xlab, "; rem: ", 
                       scales::percent(removed, accuracy = 0.1))
  }
  
  outp <- ggplot(data, aes(x = !!sym(varname_x), 
                   y = !!sym(varname_y))) +
    geom_point(alpha = alpha)
  outp <- outp + 
    geom_smooth(method = "lm", se = FALSE)  +
    facet_grid(~cond_iv_label) +
    geom_label(data = dgof, aes(label = 
                                  formatC(r, digits = 3, format = "f")), 
              x = rx[2] - 0.375*diff(rx), 
              y = ry[2] - 0.05*diff(ry), parse = TRUE, colour = "blue")
  outp +
    labs(x = new_xlab)
}

make_cor_plot2 <- function(data, x, y, filter, alpha = 0.3) {
  orig_n <- nrow(data)
  if (!missing(filter)) {
    data <- data %>% 
      dplyr::filter( !!enquo(filter) )
  }
  
  varname_x <- all.vars(enquo(x))
  varname_y <- all.vars(enquo(y))
  
  data <- data %>% 
    dplyr::filter(!is.na(!!sym(varname_x)), !is.na(!!sym(varname_y)))
  
  removed <- 1 - nrow(data)/orig_n
  
  dlm <- data %>% 
    group_by(cond_iv_label, cond_label) %>% 
    summarise(lm = list(lm(formula(expr(!!enexpr(y) ~ !!enexpr(x))))),
              r = cor({{y}}, {{x}})) 
  dlm <- dlm %>% 
    mutate(est = map(lm, ~ broom::tidy(.)), 
           gof = map(lm, ~ broom::glance(.)))
  dgof <- dlm %>% 
    select(-lm, -est) %>% 
    unnest(gof)
  rx <- range(data[[varname_x]], na.rm = TRUE)
  ry <- range(data[[varname_y]], na.rm = TRUE)
  
  new_xlab <- if (varname_x == as_label(enquo(x))) {
    varname_x
  } else {
    paste0(varname_x, " (", as_label(enquo(x)), ")")
  }
  if (!missing(filter)) {
    new_xlab <- paste0(new_xlab, "; filter: ", as_label(enquo(filter)))
  }
  if (removed > 0) {
    new_xlab <- paste0(new_xlab, "; rem: ", 
                       scales::percent(removed, accuracy = 0.1))
  }
  
  outp <- ggplot(data, aes(x = !!sym(varname_x), 
                   y = !!sym(varname_y))) +
    geom_point(alpha = alpha)
  outp <- outp + 
    geom_smooth(method = "lm", se = FALSE)  +
    facet_grid(~cond_label+cond_iv_label) +
    geom_label(data = dgof, aes(label = paste0("r = ",
                                  formatC(r, digits = 2, format = "f"))), 
              x = rx[2] - 0.375*diff(rx), 
              y = ry[2] - 0.075*diff(ry), colour = "blue")
  outp +
    labs(x = new_xlab)
}
