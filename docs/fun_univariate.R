
compare_continuous_covariate <- function(data, 
                                         covariate, 
                                         alpha = 0.3) {
  dlm <- data %>% 
    group_by(cond_x) %>% 
    summarise(lm = list(lm(formula(expr(abs_dev ~ !!enexpr(covariate)))))) 
  dlm <- dlm %>% 
    mutate(est = map(lm, ~ broom::tidy(.)), 
           gof = map(lm, ~ broom::glance(.)))
  dgof <- dlm %>% 
    select(-lm, -est) %>% 
    unnest(gof)
  
  varname <- all.vars(enquo(covariate))
  
  rx <- range(data[[varname]], na.rm = TRUE)
  
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
  
  ggplot(data, aes(x = !!sym(varname), 
                   y = abs_dev)) +
    geom_point(alpha = alpha) +
    geom_smooth(method = "lm", se = FALSE, 
                formula = formula(expr(y ~ !!expr_cov)))  +
    facet_wrap("cond_x", nrow = 1) +
    geom_label(data = dgof, aes(label = paste0("paste(italic(R) ^ 2, ' = ", 
                    substr(
                      formatC(r.squared, digits = 2, format = "f"), 
                      2, 4), "')")), 
              x = rx[2] - 0.2*diff(rx), 
              y = 0.9, parse = TRUE) +
    labs(x = new_xlab, y = ylab)
}
