
#' This file serves two purposes. Both of those concern the data files (i.e.,
#' all non-info files). First, it resaves all data files with strong
#' compression, reducing the overall size of the data. Second,  it removes data
#' attributes from all 'multiverseMPT' objects. This ensures that the data used
#' for fitting is not accidentally saved on github.

library("MPTmultiverse")

all_files <- list.files("data/", pattern = ".RData", 
                        full.names = TRUE, recursive = TRUE)
all_files <- all_files[!stringr::str_detect(all_files, "_info\\.RData$")]

for (i in seq_along(all_files)) {
  #resave <- FALSE
  objs <- load(all_files[i])
  for (j in seq_along(objs)) {
    #str(get(objs[1]), 1)
    #dplyr::glimpse(get(objs[1]))
    if (inherits(get(objs[j]), "multiverseMPT") && 
        "data" %in% names(attributes(get(objs[j])))) {
      #browser()
      #resave <- TRUE
      tmp <- get(objs[j])
      attr(tmp, "data") <- NULL
      assign(objs[j], tmp)
    }
  }
  save(list = objs, file = all_files[i], compress = "xz")
  rm(list = objs)
}
