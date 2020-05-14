
compare_continuous_covariate <- function(data, 
                                         covariate, 
                                         filter,
                                         alpha = 0.3) {
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
    group_by(cond_x) %>% 
    summarise(lm = list(lm(formula(expr(abs_dev ~ !!enexpr(covariate)))))) 
  dlm <- dlm %>% 
    mutate(est = map(lm, ~ broom::tidy(.)), 
           gof = map(lm, ~ broom::glance(.)))
  dgof <- dlm %>% 
    select(-lm, -est) %>% 
    unnest(gof)
  
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
    outp <- outp + geom_smooth(method = "gam", se = FALSE, 
                formula = y ~ s(x, bs = "cs"), colour = "red")
  }
  outp + 
    geom_smooth(method = "lm", se = FALSE, 
                formula = formula(expr(y ~ !!expr_cov)))  +
    facet_wrap("cond_x", nrow = 1) +
    geom_label(data = dgof, aes(label = paste0("paste(italic(R) ^ 2, ' = ", 
                    substr(
                      formatC(r.squared, digits = 2, format = "f"), 
                      2, 4), "')")), 
              x = rx[2] - 0.2*diff(rx), 
              y = 0.9, parse = TRUE) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = new_xlab, y = ylab)
}
