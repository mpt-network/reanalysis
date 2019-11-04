
plot_pair <- function(all_pars, cond_x, cond_y, add_ccc = TRUE) {
  pars_plot_sel <- all_pars %>% 
    filter(cond_x == {{cond_x}}, cond_y == {{cond_y}})
  
  plot_text <- pars_plot_sel %>%
    group_by(cond_x, cond_y) %>% 
    summarise(ccc = format(
      CCC(x, y, na.rm = TRUE)$rho.c$est, 
      digits = 2))
  
  plot_out <- pars_plot_sel %>% 
    ggplot(aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(alpha = 0.5) + # aes(size = trials)
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_size(range = c(1, 3.5))
    if (add_ccc) {
      plot_out <- plot_out + geom_text(data=plot_text,
                aes(x = -Inf, y = Inf, hjust = -0.04, vjust = 1.2,
                    label=ccc),
                parse = TRUE, inherit.aes=FALSE, size = 5)
    }  
  plot_out
}

get_convergence <- function(results, ..., split = TRUE) {
  if (is.null(results[["model_type"]])) {
    results$model_type <- "miss"
  }
  group_var <- enquos(...)
  
  out1 <- results %>% 
    group_split(dataset, !!! group_var) %>% 
    {suppressWarnings(map_dfr(., check_set))}
  
  out2 <- results %>% 
    group_by(dataset, !!! group_var) %>%
    group_keys()
  out <- left_join(out2, out1)
  
  if (split) {
    out %>% 
      unite("model2", !!! group_var) %>% 
      group_by(model2) %>% 
      summarise_if(is.logical, mean) %>% 
      pivot_wider(-model2) %>% 
      spread(model2, value) %>% 
      return
    
  } else {
    out %>% 
      summarise_if(is.logical, mean) %>% 
      gather() %>% 
      return
  }
  

}

check_fit <- function(df) {
  if (df$package == "TB") {
    params <- df %>%
      tidyr::unnest(.data$est_group) %>%
      dplyr::select(.data$parameter, .data$core)

    tmp_convergence <- df$convergence[[1]] %>%
      dplyr::filter(.data$Rhat > getOption("MPTmultiverse")$treebugs$Rhat_max) %>%
      dplyr::mutate(parameter = MPTmultiverse:::label_parameter(.data$parameter, params),
                    core = grepl("COREPARAMETER", .data$parameter),
                    parameter = gsub("COREPARAMETER", "", x = .data$parameter))

    if (nrow(dplyr::filter(tmp_convergence, .data$core)) > 0) {
      return(FALSE)
    }

    tmp_neff <- df$convergence[[1]] %>%
      dplyr::filter(!is.na(.data$Rhat), .data$n.eff < getOption("MPTmultiverse")$treebugs$Neff_min) %>%
      dplyr::mutate(parameter = MPTmultiverse:::label_parameter(.data$parameter, params),
                    core = grepl("COREPARAMETER", .data$parameter),
                    parameter = gsub("COREPARAMETER", "", x = .data$parameter))
    if (nrow(dplyr::filter(tmp_neff, .data$core)) > 0) {
      return(FALSE)
    }
    
  } 
  return(TRUE)
}

repair_pars <- function(df) {
  df <- df %>% 
    mutate(parameter = if_else(model == "2htsm" & 
                                 parameter == "d_3", "d_2", 
                               parameter)) %>% 
    mutate(parameter = factor(paste0(model, ":", parameter)))
  df
}

get_rel_par_weight <- function(parameter, est, model_exp, rel_tree, 
                               orig_condition) {
  if (is.null(rel_tree[[1]])) return(rep(NA_real_, length(parameter)))
  
  ## check all trees are equal: 
  stopifnot(map_lgl(model_exp, ~all.equal(., model_exp[[1]])))
  
  ## calculate predictions per branch
  tmpe <- new.env()
  for (i in seq_along(parameter)) {
    assign(parameter[i], est[i], tmpe)
  }
  try(branch_prob <- map_depth(model_exp[[1]], 2, eval, envir = tmpe))
  if (!exists("branch_prob")) {
    #browser()
    return(rep(NA_real_, length(parameter)))
  }
  
  ## check if probabilities on trees sum to 1
  if (isTRUE(any(abs(map_dbl(branch_prob, ~sum(unlist(.))) - 1) > 0.01)))
    browser()
  
  rel_weight <- vector("numeric", length(parameter))
  pars_per_branch <- map_depth(model_exp[[1]], 2, all.vars)
  for (i in seq_along(rel_weight)) {
    par_in_branch <- map_depth(pars_per_branch, 2, ~parameter[[i]] %in% .) 
    par_info_tree <- map2_dbl(branch_prob, par_in_branch, 
                              ~sum(unlist(.x)[unlist(.y)]))
    rel_weight[i] <- rel_tree[[1]] %>% 
      filter(condition == orig_condition[1]) %>% 
      select(-condition) %>% 
      unlist %>% 
      {. * par_info_tree} %>% 
      sum
    
  }
  rel_weight
  
}

get_weight_n <- function(parameter, est, model_exp, data_tree, 
                               orig_condition) {
  if (is.null(data_tree[[1]])) return(rep(NA_real_, length(parameter)))
  
  ## check all trees are equal: 
  stopifnot(map_lgl(model_exp, ~all.equal(., model_exp[[1]])))
  
  ## calculate predictions per branch
  tmpe <- new.env()
  for (i in seq_along(parameter)) {
    assign(parameter[i], est[i], tmpe)
  }
  try(branch_prob <- map_depth(model_exp[[1]], 2, eval, envir = tmpe))
  if (!exists("branch_prob")) return(rep(NA_real_, length(parameter)))
  
  ## check if probabilities on trees sum to 1
  if (isTRUE(any(abs(map_dbl(branch_prob, ~sum(unlist(.))) - 1) > 0.01)))
    browser()
  
  rel_weight <- vector("numeric", length(parameter))
  pars_per_branch <- map_depth(model_exp[[1]], 2, all.vars)
  for (i in seq_along(rel_weight)) {
    par_in_branch <- map_depth(pars_per_branch, 2, ~parameter[[i]] %in% .) 
    par_info_tree <- map2_dbl(branch_prob, par_in_branch, 
                              ~sum(unlist(.x)[unlist(.y)]))
    rel_weight[i] <- data_tree[[1]] %>% 
      group_by(condition) %>% 
      select(-id) %>% 
      summarise_all(sum) %>% 
      filter(condition == orig_condition[1]) %>% 
      select(-condition) %>% 
      unlist %>% 
      {. * par_info_tree} %>% 
      sum
    
  }
  rel_weight
  
}

fit_polynomials <- function(formula, data, ..., max_degree = 8) {
  out <- vector("list", max_degree)
  for (i in seq_len(max_degree)) {
    tmpf <- do.call("substitute", list(formula, list(DEG=i)))
    out[[i]] <- lmer(formula = tmpf, data = data, ...)
    out[[i]]@call$formula <- tmpf
    out[[i]]@call$data <- substitute(data)
  }
  out
}
