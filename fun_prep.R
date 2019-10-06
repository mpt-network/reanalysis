
load_combine <- function(files) {
  tmpe <- map(files, ~new.env())
  ## load file names in separate environment
  obj <- map2_chr(files, tmpe, ~load(.x, envir = .y))
  ## check if dataset column is empty and if so, set from file name
  dnames <- str_remove_all(string = files, pattern = ".+/") %>% 
    str_remove_all("\\.RData$") %>% 
    str_remove_all("\\.rda$")
  try({
  for (i in seq_along(obj)) {
    if (tmpe[[i]][[obj[i]]]$dataset[1] == "") {
      cat(dnames[i], " ", obj[i], "\n")
      tmpe[[i]][[obj[i]]]$dataset <- dnames[i]
    }
  }
  }, silent = TRUE)
  if (!is.data.frame(tmpe[[1]][[obj[1]]])) {
    results <- map2_dfr(obj, tmpe, ~{
      out <- get(.x, envir = .y)
      out <- map(out, prepare_single_result)
      bind_rows(out)
    })
  } else {
    results <- map2_dfr(obj, tmpe, ~{
      out <- get(.x, envir = .y)
      prepare_single_result(out)
    })
  }
  results$orig_model <- results$model
  sanitize_vars(results)
}

load_combine_info <- function(files) {
  tmpe <- map(files, ~new.env())
  ## load file names in separate environment
  # obj <- map2(files, tmpe, ~list(load(.x, envir = .y))) %>% 
  #   unlist
  obj <- map2_chr(files, tmpe, ~load(.x, envir = .y)) 
  
  ## check if dataset column is empty and if so, set from file name
  # dnames <- str_remove_all(string = files, pattern = ".+/") %>% 
  #   str_remove_all("\\.RData$") %>% 
  #   str_remove_all("\\.rda$") %>% 
  #   str_remove_all("_info$")
  # names <- str_split(dnames, "\\.csv_", simplify = TRUE)
  # colnames(names) <- c("data", "orig_model")
  # names <- as_tibble(names)
  
  out <- map2_dfr(obj, tmpe, make_info_df)
  # browser()
  out
}

make_info_df <- function(object, env) {
  out <- get(object, envir = env)
  colnames(out$data_tree)[1] <- "id"
  colnames(out$data_tree)[2] <- "condition"
  n_tree <- out$data_tree %>% 
    group_by(condition) %>% 
    select(-1) %>% 
    summarise_all(sum) %>% 
    ungroup
  rel_tree <- n_tree
  rel_tree[,-1] <- rel_tree[,-1]/rowSums(rel_tree[,-1])
  # str(out, 1)
  df <- tibble(
    dataset = out$dataset,
    orig_model = out$model,
    n_tree = list(n_tree),
    rel_tree = list(rel_tree),
    model_df = list(out$model_eqn),
    model_exp = list(out$model_branches),
    data_tree = list(out$data_tree)
  )
  return(df)
}

prepare_single_result <- function(result) {
  out <- result
  if (!is.null(attr(out, "data", exact = TRUE))) {
    trials <- tibble(
      id = attr(out, "data")[[1]],
      condition = attr(out, "data")[[2]],
      trials = rowSums(attr(out, "data")[,-(1:2)])
    ) 
    out$trials <- list(trials) 
  } else {
    out$trials <- list(NULL)
  }
  
  for (i in seq_len(nrow(out))) {
    out[i,]$est_group[[1]]$orig_condition <- out[i,]$est_group[[1]]$condition
    out[i,]$est_group[[1]]$orig_parameter <- out[i,]$est_group[[1]]$parameter
    if (nrow(out[i,]$est_indiv[[1]]) > 0) {
      out[i,]$est_indiv[[1]]$orig_condition <- out[i,]$est_indiv[[1]]$condition
      out[i,]$est_indiv[[1]]$orig_parameter <- out[i,]$est_indiv[[1]]$parameter
    }
  }
  out
}

check_core_pars <- function(results) {
  all_core <- results %>%
    unnest(est_group) %>%
    group_by(model) %>%
    group_split() %>%
    map(~unique(.$parameter[.$core])) %>% 
    unlist %>% 
    unique() %>% 
    na.omit
  all_non_core <- results %>%
    unnest(est_group) %>%
    group_by(model) %>%
    group_split() %>%
    map(~unique(.$parameter[!.$core])) %>% 
    unlist %>% 
    unique() %>% 
    na.omit()
  na_parameters <- results %>%
    unnest(est_group) %>%
    group_by(model) %>%
    group_split() %>%
    map(~unique(.$parameter[is.na(.$core)])) %>% 
    unlist %>% 
    unique() %>% 
    na.omit()
  if (length(all_core) > 0)
    cat("Core parameters: \n", all_core, "\n")
  if (length(all_non_core) > 0)
    cat("Non-core parameters: \n", all_non_core, "\n")
  if (length(na_parameters) > 0)
    cat("NA parameters: \n", na_parameters, "\n")
  if (any(all_core %in% all_non_core)) {
    cat("Parameters core and non-core:\n", 
        all_core[all_core %in% all_non_core], "\n")
  }
}

check_beta <- function(results) {
  split_res <- results %>% 
    group_by(dataset, model)
  cat("N datasets:", length(group_split(split_res)), "\n")
  
  beta_ns <- split_res %>%
    summarise(n_beta = sum(method == "beta"),
              n_betacpp = sum(method == "beta++"))
  
  try(prior_beta <- split_res %>% 
    filter(method == "beta") %>% 
    select(model:method, options) %>% 
    unnest(options) %>% 
    select(dataset, model, prior.beta))
  if (!exists("prior_beta")) {
    split_res$options <- map(split_res$options, ~{
      .[["n.chains"]] <- as.numeric(.[["n.chains"]])
      .})
    prior_beta <- split_res %>% 
      filter(method == "beta") %>% 
      select(model:method, options) %>% 
      unnest(options) %>% 
      select(dataset, model, prior.beta)
  }
  
  out <- left_join(beta_ns, prior_beta, by = c("dataset", "model")) 
  out %>% 
    print(n = 1e6)
  invisible(out)
}

check_double_fit <- function(results) {
  if (is.null(results[["model_type"]])) {
    results$model_type <- "miss"
  }
  double <- results %>% 
    group_by(dataset, model, model_type, inter) %>% 
    tally() %>% 
    filter(n != 1)
  if (nrow(double) > 0) {
    cat("Encountered cases with more than one fit:\n")
    print(double, n = 1e5)
  }
  
  check_double <- vector("list", nrow(double))
  for (i in seq_along(check_double)) {
    tmp_res <- left_join(double[i,], results, 
              by = c("dataset", "model", "model_type", "inter")) %>% 
      ungroup
    
    ## if bayes
    if (tmp_res[1,]$package == "TB") {
      params <- tmp_res[1,] %>%
        tidyr::unnest(.data$est_group) %>%
        dplyr::select(.data$parameter, .data$core)
      
      check_tmp <- tibble(
        n = seq_len(nrow(tmp_res)),
        rhat = TRUE,
        neff = TRUE
      )
      for (j in seq_len(nrow(tmp_res))) {
        tmp_convergence <- tmp_res[j,]$convergence[[1]] %>%
          dplyr::filter(.data$Rhat > getOption("MPTmultiverse")$treebugs$Rhat_max) %>%
          dplyr::mutate(parameter = 
                          MPTmultiverse:::label_parameter(.data$parameter, params),
                        core = grepl("COREPARAMETER", .data$parameter),
                        parameter = gsub("COREPARAMETER", "", x = .data$parameter))
        
        if (nrow(dplyr::filter(tmp_convergence, .data$core)) > 0) {
          check_tmp[j , "rhat"  ] <- FALSE
        }
        
        tmp_neff <- tmp_res[j,]$convergence[[1]] %>%
          dplyr::filter(!is.na(.data$Rhat), .data$n.eff < 
                          getOption("MPTmultiverse")$treebugs$Neff_min) %>%
          dplyr::mutate(parameter = 
                          MPTmultiverse:::label_parameter(.data$parameter, params),
                        core = grepl("COREPARAMETER", .data$parameter),
                        parameter = gsub("COREPARAMETER", "", x = .data$parameter))
        if (nrow(dplyr::filter(tmp_neff, .data$core)) > 0) {
          check_tmp[j , "neff"  ] <- FALSE
        }
      }
      check_double[[i]] <- tmp_res[which.max(rowSums(check_tmp[,-1]))[1],]
    } else {
      check_double[[i]] <- tmp_res[nrow(tmp_res),]
    }
  }
  out <- anti_join(results, double,
              by = c("dataset", "model", "model_type", "inter")) %>% 
    bind_rows(check_double)
  out$n <- NULL
  out
}



get_convergence <- function(results) {
  if (is.null(results[["model_type"]])) {
    results$model_type <- "miss"
  }
  
  out1 <- results %>% 
    group_split(dataset, model, model_type) %>% 
    {suppressWarnings(map_dfr(., check_set))}
  
  out2 <- results %>% 
    group_by(dataset, model, model_type) %>%
    group_keys()
  out <- left_join(out2, out1, by = c("dataset", "model"))
  
  out %>% 
    unite("model2", model, model_type) %>% 
    group_by(model2) %>% 
    summarise_if(is.logical, mean) %>% 
    pivot_longer(-model2) %>% 
    spread(model2, value)

}

# get_info <- function(list) {
#   cat("Models: \n")
#   map_chr(list, ~unique(.$model)) %>% 
#     unique %>% 
#     sort %>% 
#     print
#   
#   cat("Data: \n")
#   map_chr(list, ~unique(.$dataset)) %>% 
#     unique %>% 
#     sort %>% 
#     print
#   
#   cat("Method: \n")
#   map(list, ~unique(unite(select(., pooling, package, method), col = "meth")$meth)) %>%
#     unlist %>% 
#     unique %>% 
#     sort %>% 
#     print
# }

get_info_df <- function(df) {
  cat("Models: \n")
  df$model %>% 
    unique %>% 
    sort %>% 
    print

  cat("Data: \n")
  df$dataset %>% 
    unique %>% 
    sort %>% 
    print
  
  cat("Method: \n")
  df %>% 
    select(pooling, package, method) %>% 
    unite(col = "meth") %>% 
    {unique(.$meth)} %>% 
    sort %>% 
    print
}


check_set <- function(results) {

  if (!all(results$model[1] == results$model))
    stop("All results need to use same model")
  if (!all(results$dataset[1] == results$dataset))
    stop("All results need to use same dataset")

  # expected <- structure(list(
  #   pooling = c("no", "no", "no", "complete", "no", "complete", "partial",
  #               "partial", "partial", "partial", "partial"),
  #   package = c("MPTinR", "MPTinR", "MPTinR", "MPTinR", "TreeBUGS", "TreeBUGS",
  #               "TreeBUGS", "TreeBUGS", "TreeBUGS", "TreeBUGS", "HMMTreeR"),
  #   method = c("NPB/MLE", "PB/MLE", "asymptotic", "asymptotic", "simple",
  #              "simple", "trait", "trait_uncorrelated", "beta", "betacpp", 
  #              "latent_class")),
  #   .Names = c("pooling", "package", "method"),
  #   class = c("tbl_df", "tbl", "data.frame"
  #   ), row.names = c(NA, -11L))
  expected <- structure(list(
    pooling = c("No", "No", "No", "Comp", "No", "Comp", "PP",
                "PP", "PP", "PP", "PP"),
    package = c("MR", "MR", "MR", "MR", "TB", "TB",
                "TB", "TB", "TB", "TB", "LC"),
    method = c("PB", "NPB", "MLE", "MLE", 
               "ss", "ss",
               "trait", "trait_u", "beta", "beta++", "lc")),
    .Names = c("pooling", "package", "method"),
    class = c("tbl_df", "tbl", "data.frame"
    ), row.names = c(NA, -11L))
  missing <- dplyr::anti_join(expected, results[, 3:5], by = c("pooling", "package", "method"))
  miss <- tidyr::unite(missing, "meth")

  out <- tidyr::unite(expected, "meth")
  out$problem <- TRUE

  for (i in seq_len(nrow(miss))) {
    out[ out$meth == miss$meth[i] , "problem"  ] <- FALSE
  }

  ### TreeBUGS
  res_tree <- results %>%
    dplyr::filter(.data$package == "TB") %>%
    dplyr::select(!!c("model", "dataset", "pooling", "package", "method", "convergence", "est_group"))

  for (i in seq_len(nrow(res_tree))) {

    cur_meth <- tidyr::unite(res_tree[i,3:5], "meth")$meth

    params <- res_tree[i,] %>%
      tidyr::unnest(.data$est_group) %>%
      dplyr::select(.data$parameter, .data$core)

    tmp_convergence <- res_tree[i, ]$convergence[[1]] %>%
      dplyr::filter(.data$Rhat > getOption("MPTmultiverse")$treebugs$Rhat_max) %>%
      dplyr::mutate(parameter = MPTmultiverse:::label_parameter(.data$parameter, params),
                    core = grepl("COREPARAMETER", .data$parameter),
                    parameter = gsub("COREPARAMETER", "", x = .data$parameter))

    if (nrow(dplyr::filter(tmp_convergence, .data$core)) > 0) {
      out[ out$meth == cur_meth , "problem"  ] <- FALSE
    }

    tmp_neff <- res_tree[i,]$convergence[[1]] %>%
      dplyr::filter(!is.na(.data$Rhat), .data$n.eff < getOption("MPTmultiverse")$treebugs$Neff_min) %>%
      dplyr::mutate(parameter = MPTmultiverse:::label_parameter(.data$parameter, params),
                    core = grepl("COREPARAMETER", .data$parameter),
                    parameter = gsub("COREPARAMETER", "", x = .data$parameter))
    if (nrow(dplyr::filter(tmp_neff, .data$core)) > 0) {
      out[ out$meth == cur_meth , "problem"  ] <- FALSE
    }
  }
  dplyr::bind_cols(
    results[1, 1:2]
    ,
    tidyr::spread(out, "meth", "problem")
  )
}

sanitize_vars <- function(results) {
  
  results$pooling <- factor(results$pooling, 
                            levels = c("no", "complete",  "partial"), 
                            labels = c("No", "Comp", "PP"))
  
  results$package <- factor(results$package, 
                            levels = c("MPTinR", "TreeBUGS", "HMMTreeR"), 
                            labels = c("MR", "TB", "LC"))
  results$method <- factor(results$method, 
                           levels = c("asymptotic", "PB/MLE", "NPB/MLE", 
                                      "simple", 
                                      "trait", "trait_uncorrelated",
                                      "betacpp", "beta", "latent_class"), 
                           labels = c("MLE", "PB", "NPB", "ss", 
                                      "trait", "trait_u",
                                      "beta++", "beta", "lc"))
  
  results$inter <- with(results, interaction(method, pooling, package, 
                                             drop = TRUE, sep = " "))
  results$inter <- fct_recode(results$inter,
                              "Comp MLE" = "MLE Comp MR", 
                              "Comp Bayes" = "ss Comp TB",
                              "No PB" = "PB No MR", 
                              "No NPB" = "NPB No MR",
                              "No asy" = "MLE No MR", 
                              "No Bayes" = "ss No TB",
                              "Trait PP" = "trait PP TB", 
                              "Trait_u PP" = "trait_u PP TB",
                              "Beta PP" = "beta PP TB",
                              "Beta++ PP" = "beta++ PP TB",
                              "LC PP" = "lc PP LC")
  results$inter <- fct_relevel(results$inter, 
                               c("Comp MLE", "Comp Bayes", 
                                 "No asy", "No PB", "No NPB", "No Bayes",  
                                 "LC PP",
                                 "Beta PP", "Beta++ PP", 
                                 "Trait_u PP", "Trait PP"))
  results
}


