
make_categorical_table <- function(data, 
                                   covariate, 
                                   fun = mean,
                                   ...,
                                   filter, 
                                   add_prop = TRUE) {
  orig_n <- nrow(data)
  if (!missing(filter)) {
    data <- data %>% 
      dplyr::filter( !!enquo(filter) )
  }
  removed <- 1 - nrow(data)/orig_n
  
  tmp_out <- data %>% 
    group_by(!!enquo(covariate), cond_x) %>% 
    summarise(tmp = fun(abs_dev, ...)) %>%
    pivot_wider(names_from = cond_x, values_from = tmp)
  
  if (add_prop) {
    tmp_prop <-  data %>% 
      group_by(!!enquo(covariate)) %>% 
      summarise(prop = n()/nrow(.))
    tmp_out <- left_join(tmp_prop, tmp_out)
  }
  tmp_out %>% 
    print(n = Inf)
}

compare_categorical_covariate <- function(data, 
                                          covariate, 
                                          filter, 
                                          plot = TRUE,
                                          max_level_plot = 9,
                                          rm_levels_less = 0.04,
                                          alpha = 0.3) {
  
  if (!missing(filter)) {
    prep_data_model <- fit_model_covariate(data = data, 
                                           covariate = enquo(covariate), 
                                           filter = enquo(filter))
  } else {
    prep_data_model <- fit_model_covariate(data = data,
                                           covariate = enquo(covariate))
  }
  
  props <- prep_data_model$data %>% 
    group_by(!!enquo(covariate)) %>% 
    tally() %>% 
    mutate(prop = n / sum(n))
  
  if (nrow(props) > max_level_plot) {
    levels_use <- props[[1]][props$prop >= rm_levels_less]
  } else {
    levels_use <- props[[1]]
  }
  
  if (plot && length(levels_use) > 0 && length(levels_use) <= max_level_plot) {
    new_xlab <- prep_data_model$varname
    if (length(levels_use) < nrow(props)) {
      new_xlab <- paste0(new_xlab, "; rem: ", 
                         nrow(props) - length(levels_use), " levels")
    }
    pout <- data %>% 
      dplyr::filter(!!sym(prep_data_model$varname) %in% levels_use) %>% 
      ggplot(aes(x = !!sym(prep_data_model$varname), 
                   y = abs_dev)) +
      geom_violin() + 
      stat_summary(fun.data = mean_se) +
      facet_wrap("cond_x", nrow = 1) +
    geom_label(
      data = prep_data_model$dgof, 
               aes(label = paste0("paste(italic(R) ^ 2, ' = ", 
                    substr(
                      formatC(r.squared, digits = 2, format = "f"), 
                      2, 4), "')")), 
              x = max(1.5, length(levels_use)-1.5), 
              y = 0.9, parse = TRUE) +
    coord_cartesian(ylim = c(0, 1)) +
      labs(x = new_xlab, y = ylab) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(pout)
  }
  
  df_out <- prep_data_model$dgof %>% 
    select(cond_x, r.squared) %>% 
    mutate(covariate = prep_data_model$varname, 
           nlevels = nrow(props)) %>% 
    pivot_wider(id_cols = c(covariate, nlevels), 
                names_from = cond_x, 
                values_from = r.squared)
  return(df_out)
  
}

fit_model_covariate <- function(data, 
                                covariate, 
                                filter) {
  orig_n <- nrow(data)
  if (!missing(filter)) {
    data <- data %>% 
      dplyr::filter( !!enquo(filter) )
  }
  
  varname <- all.vars(covariate)
  
  data <- data %>% 
    dplyr::filter(!is.na(!!sym(varname)))
  
  removed <- 1 - nrow(data)/orig_n

  dlm <- data %>% 
    group_by(cond_x) %>% 
    nest() %>% 
    mutate(lm = map(data, 
                    ~lm(formula(expr(abs_dev ~ !!ensym(covariate))), .)))
  dlm <- dlm %>% 
    mutate(est = map(lm, ~ broom::tidy(.)), 
           gof = map(lm, ~ broom::glance(.)))
  dgof <- dlm %>% 
    select(-lm, -est) %>% 
    unnest(gof)
  return(list(
    data = data,
    dlm = dlm, 
    dgof = dgof, 
    varname = varname, 
    removed = removed
  ))
}

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
                formula = y ~ s(x, bs = "ts"), colour = "red")
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
