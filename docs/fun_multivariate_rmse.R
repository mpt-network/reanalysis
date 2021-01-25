
lm_scale <- function(formula, ...) {
  mmtmp <- model.frame(formula(formula), data = parent.frame())
  #head(mmtmp)
  for (i in 1:ncol(mmtmp)) {
    mmtmp[[i]] <- unlist(scale(mmtmp[[i]]))
  }
  mmtmp <- as.matrix(mmtmp)
  
  lm_out <- lm(mmtmp[,1] ~ mmtmp[,-1])
  attr(lm_out$coefficients, "names") <- 
    str_remove(attr(lm_out$coefficients, "names"), "mmtmp\\[, -1\\]")
  lm_out
}

create_filter_from_formula <- function(formula) {
  mm <- model.frame(formula(formula), data = parent.frame(), 
                    na.action = NULL)
  unname(!apply(mm, 1, function(x) any(is.na(x))))
}

produce_interaction_matrix <- function(inter_ind, length) {
  mat <- matrix(FALSE, nrow = length, ncol = length)
  for (i in inter_ind) {
    mat[ as.numeric(substr(i, 2, 2)), as.numeric(substr(i, 4, 4))  ] <- TRUE
  }
  mat
}

fit_ramp <- function(formula, center = FALSE) {
  mmtmp <- model.frame(formula(formula), data = parent.frame())
  #head(mmtmp)
  #browser()
  colnames(mmtmp) <- make.names(colnames(mmtmp))
  if (center) {
    for (i in seq_len(ncol(mmtmp)-1)) {
      if (colnames(mmtmp)[i+1] == "y_c") next
      mmtmp[,(i+1)] <- mmtmp[,(i+1)] - mean(mmtmp[,(i+1)]) 
    }
  }
  rout <- tryCatch(RAMP(X = mmtmp[,-1], y = mmtmp[, 1]), error = function(e) NA)
  rout
}

get_rmse_ramp <- function(model, data) {
  data$predict <- predict(model)
  sqrt(mean((data$abs_dev - data$predict)^2))
}

# best_two <- function(global.model, fixed) {
#   browser()
#   dd <- dredge(global.model = global.model, fixed = fixed)
# }

best_x <- function(data, formula, fixed, x = 2) {
  mod <- lm(formula(formula), data = data)
  dd <- dredge(global.model = mod, fixed = fixed, m.lim = c(x,x), 
               rank = sigma)
  mod_out <- eval(getCall(dd, 1))
  mod_out
}
