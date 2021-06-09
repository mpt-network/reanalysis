
compare_continuous_covariate <- function(data, 
                                         covariate, 
                                         ...,
                                         filter,
                                         ylab,
                                         alpha = 0.3, 
                                         scales = "fixed") {
  orig_n <- nrow(data)
  if (!missing(filter)) {
    data <- data %>% 
      dplyr::filter( !!enquo(filter) )
  }
  
  varname <- all.vars(enquo(covariate))
  
  data <- data %>% 
    dplyr::filter(!is.na(!!sym(varname)))
  
  removed <- 1 - nrow(data)/orig_n

  dlm <- data %>% 
    group_by(...) %>% 
    summarise(lm = list(lm(formula(expr(abs_dev ~ !!enexpr(covariate))))),
              xpos_label1 = max(!!sym(varname)) - 0.4*diff(range(!!sym(varname))),
              xpos_label2 = max(!!sym(varname)) - 0.12*diff(range(!!sym(varname)))) 
  dlm <- dlm %>% 
    mutate(est = map(lm, ~ broom::tidy(.)), 
           gof = map(lm, ~ broom::glance(.)))
  dgof <- dlm %>% 
    select(-lm, -est) %>% 
    unnest(gof)
  if (scales == "fixed") {
    dgof$xpos_label1 <- max(data[[varname]]) - 0.4*diff(range(data[[varname]]))
    dgof$xpos_label2 <- max(data[[varname]]) - 0.12*diff(range(data[[varname]]))
  }
  
  #rx <- range(data[[varname]], na.rm = TRUE)
  
  
  # rx <- data %>% 
  #   summarise(range = list(range(!!enquo(covariate)))) %>% 
  #   {.$range[[1]]}
  # data <- data %>% 
  #   mutate(new_x = as.vector(!!enquo(covariate)))
  
  varname <- all.vars(enquo(covariate))
  enquo_cov <- enquo(covariate)
  expr_cov <- rlang::get_expr(enquo_cov)
  if (is.call(expr_cov)) {
    expr_cov[[which(map_chr(rlang::call_args(expr_cov), as_label) == varname) + 
                1]] <- quote(x)
  } else {
    expr_cov <- quote(x)
  }
  
  new_xlab <- if (varname == as_label(enquo(covariate))) {
    varname
  } else {
    paste0(varname, " (", as_label(enquo(covariate)), ")")
  }
  if (!missing(filter)) {
    new_xlab <- paste0(new_xlab, "; filter: ", as_label(enquo(filter)))
  }
  if (removed > 0) {
    new_xlab <- paste0(new_xlab, "; rem: ", 
                       scales::percent(removed, accuracy = 0.1))
  }
  
  outp <- ggplot(data, aes(x = !!sym(varname), 
                   y = abs_dev)) +
    geom_point(alpha = alpha)
  if (INCLUDE_GAM) {
    gamr2 <- function(gam){
      1-((sum(residuals(gam)^2))/
                 (sum((gam$y - mean(gam$y))^2)))
    }
    gam_sigma <- function(gam) {
      sqrt(mean(residuals(gam)^2))
    }
    dgam <- data %>% 
      group_by(...) %>% 
      summarise(lm = list(mgcv::gam(formula(expr(
        abs_dev ~ s(!!enexpr(covariate), bs = "ts")))))) 
    dgam <- dgam %>% 
      mutate(r.squared = map_dbl(lm, gamr2),
             sigma = map_dbl(lm, gam_sigma),
             est = map(lm, ~ broom::tidy(.)), 
             gof = map(lm, ~ broom::glance(.)))
    dgam <- left_join(dgam, select(dgof, cond_label, cond_iv_label, xpos_label2))
    outp <- outp + geom_smooth(method = "gam", se = FALSE, 
                formula = y ~ s(x, bs = "ts"), colour = "red")
  }
  outp <- outp + 
    geom_smooth(method = "lm", se = FALSE, 
                formula = formula(expr(y ~ !!expr_cov)))  +
    facet_wrap(vars(...), nrow = 1, scales = scales) +
    geom_label(data = dgof, aes(
      label = paste0("'", #paste0("paste(bar(Delta), ' = ", 
                    substr(
                      formatC(sigma, digits = 3, format = "f"), 
                      2, 5), "'"), 
      x = xpos_label1), #"')")), 
              #x = rx[2] - 0.4*diff(rx), 
              y = 0.94, parse = TRUE, colour = "blue")
  if (INCLUDE_GAM) {
    outp <- outp + 
      geom_label(data = dgam, aes(label = paste0(#"paste(italic(R) ^ 2, ' = ", 
                                                 substr(
                      formatC(sigma, digits = 3, format = "f"), 
                      2, 5)), 
                      x = xpos_label2),#"')")), 
              #x = rx[2] - 0.12*diff(rx), 
              y = 0.94, parse = FALSE, colour = "red")
  }
  
  outp +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = new_xlab, y = ylab)
}

create_filter_from_formula <- function(formula) {
  mm <- model.frame(formula(formula), data = parent.frame(), 
                    na.action = NULL)
  unname(!apply(mm, 1, function(x) any(is.na(x))))
}


make_cut_labels <- function(labels) {
  nums <- map(str_extract_all(labels, "0\\.?\\d*"), as.numeric)
  substr(as.character(map_dbl(nums, 1) + map_dbl(nums, diff)/2), 2, 10)
}
