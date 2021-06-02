
harmonize_parname <- function(model, model2, parameter) {
  case_when(
    model2 == "2htsm_4" ~ substr(parameter, 1, 1),
    
    model2 %in% c("2htsm_5d", "2htsm_6e") & str_detect(parameter, "_") ~ 
      paste0(substr(parameter, 1, 1), "_", 
             substr(parameter, nchar(parameter), nchar(parameter))),
    model2 %in% c("2htsm_5d", "2htsm_6e") & !str_detect(parameter, "_") ~ 
      substr(parameter, 1, 1),
    
    model == "pc" ~ substr(parameter, 1, 1),
    
    model == "pd" ~ rename_pd_pars(parameter),
    
    model == "pm" ~ str_remove(parameter, "[[:lower:]]+"),
    
    TRUE ~ parameter
  )
}

rename_pd_pars <- function(parameter) {
  new_par_names <- c(
    a1 = "A" #1
    , a2 = "A" #2
    , a3 = "A" #3
    , c1 = "C" #1
    , c2 = "C" #2
    , c3 = "C" #3
    , c = "C" #1
    , u = "A" #1
    , a = "A" #1
    , H_50 = "A"
    , H_75 = "A"
    , R_50 = "C"
    , R_75 = "C"
    , h1 = "A_alt"
    , h2 = "A_alt"
    , h3 = "A_alt"
    , wi = "C_Exclusion"
    , wc = "C_Inclusion"
    , A = "A"
    , C = "C"
    , ab1 = "A_alt" # ws1
    , ab2 = "A_alt" # ws2
    , ab3 = "A_alt" # ws3
    , aw1 = "A_alt" # ws4
    , aw2 = "A_alt" # ws5
    , aw3 = "A_alt" # ws6
    , bg = "C_Inclusion"
    , bt = "C_Exclusion"
    , wg = "C_Exclusion"
    , wt = "C_Inclusion"
    , xx = "xx"
  )
  parameter[!(parameter %in% names(new_par_names))] <- "xx"
  new_par_names[parameter]
}

get_red_corr <- function(model, model2, parameter1, parameter2, parameter_o, 
                         variable, operation, neutral = NA_real_) {
  stopifnot(all(parameter1[1] == parameter1))
  stopifnot(all(model2[1] == model2))
  stopifnot(all(model[1] == model))
  
  if (model[1] == "2htsm") {
    if (parameter1[1] %in% c("d", "d_1", "d_2")) {
      variable_use <- variable[parameter2 != "b"]
    } else if (parameter1[1] == "b") {
      variable_use <- variable[!(parameter2 %in% c("d", "d_1", "d_2"))]
    } else if (parameter1[1] %in%  c("D", "D_1", "D_2","g")) {
      variable_use <- variable
    }
  }
  
  if (model[1] == "c2ht") {
    if (parameter1[1] %in% c("Dn")) {
      variable_use <- variable[!(parameter2 %in% c("Do"))]
    } else if (parameter1[1] == "Do") {
      variable_use <- variable[!(parameter2 %in% c("Dn"))]
    } else if (parameter1[1] %in%  c("g")) {
      variable_use <- variable[!(parameter2 %in% 
                                   c("r_1", "r_2", "r_6", "r_3", "r_8"))]
    }
  }
  
  if (model[1] == "pc") {
    if (parameter1[1] %in% c("r")) {
      variable_use <- variable[parameter2 != "u"]
    } else if (parameter1[1] == "u") {
      variable_use <- variable[!(parameter2 %in% c("r"))]
    } else if (parameter1[1] %in%  c("c")) {
      variable_use <- variable
    }
  }
  
  if (model[1] == "pd") {
    if (parameter1[1] %in% c("C_Inclusion")) {
      variable_use <- variable[parameter2 != "C_Exclusion"]
    } else if (parameter1[1] %in% c("C_Exclusion")) {
      variable_use <- variable[parameter2 != "C_Inclusion"]
    } else if (parameter1[1] %in%  c("A", "C", "A_alt")) {
      variable_use <- variable
    }
  }
  
  if (model[1] == "pm") {
    if (parameter1[1] %in% c("C1")) {
      variable_use <- variable[!(parameter2 %in% c("C2"))]
    } else if (parameter1[1] == "C2") {
      variable_use <- variable[!(parameter2 %in% c("C1"))]
    } else if (parameter1[1] %in%  c("M", "P")) {
      variable_use <- variable
    }
  }
  
  if (model[1] == "hb") {
    if (parameter1[1] %in% c("rc")) {
      variable_use <- variable[!(parameter2 %in% 
                                   c("le", "re","gg3", "gl3", "b", "c"))]
    } else if (parameter1[1] == "re") {
      variable_use <- variable[!(parameter2 %in% c("lc", "rc"))]
    } else if (parameter1[1] == "b") {
      variable_use <- variable[!(parameter2 %in% c("lc", "rc"))]
    } else if (parameter1[1] == "c") {
      variable_use <- variable[!(parameter2 %in% 
                                   c("lc", "rc", "gl1", "gl2", "gg1", "gg2"))]
    }
  }
  
  if (model[1] == "rm") {
    if (parameter1[1] %in% c("g")) {
      return(neutral)
    } else if (parameter1[1] %in% c("a", "r", "b")) {
      variable_use <- variable[!(parameter2 %in% c("g"))]
    } 
  }
  
  if (model[1] == "real") {
    if (parameter1[1] %in% c("A1", "L1")) {
      variable_use <- variable[(parameter2 %in% 
                                  c("A1" ,"Re", "attReC", "attReT", "L1", "attL"))]
    } else if (parameter1[1] %in% c("A2", "L2")) {
      variable_use <- variable[(parameter2 %in% 
                                  c("A2", "Re", "attReC", "attReT", "L2", "attL"))]
    } else if (parameter1[1] %in% c("L3")) {
      variable_use <- variable[(parameter2 %in% 
                                  c("Re", "attReT", "attL"))]
    } else if (parameter1[1] %in% c("L4")) {
      variable_use <- variable[(parameter2 %in% 
                                  c("Re", "attReT", "attL"))]
    } else if (parameter1[1] %in%  c("Re")) {
      variable_use <- variable
    }
  }
  
  if (model[1] == "quad") {
    if (parameter1[1] %in% c("ACbb1")) {
      variable_use <- variable[!(parameter2 %in% c("ACwg1"))]
    } else if (parameter1[1] %in% c("ACwg1")) {
      variable_use <- variable[!(parameter2 %in% c("ACbb1"))]
    } else if (parameter1[1] %in% c("G1")) {
      variable_use <- variable[!(parameter2 %in% c("OB1"))]
    } else if (parameter1[1] %in% c("OB1")) {
      variable_use <- variable[!(parameter2 %in% c("G1"))]
    } else if (parameter1[1] %in%  c("D1")) {
      variable_use <- variable
    }
  }
  
  if (exists("variable_use") && length(variable_use) == 0) {
    browser()
    print(parameter_o[1])
  }
  
  tryCatch(operation(variable_use), error = function(e) NA_real_)
}
  
