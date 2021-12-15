#' Deprecated. See saveLoad__saveCompactLooks
#'
#' @export
saveCompactLooks = function(...){ saveLoad__saveCompactLooks(...) }

#' Save all the leave-one-out knockoffs to disk, e.g. for reading into Python.
#'
#' Save them in their low-rank representation to avoid wasting space.
#'
#' @export
saveLoad__saveCompactLooks = function(looks, savepath, filetype = c("csv", "h5")){
  dir.create(savepath, recursive = T, showWarnings = F)
  filetype = match.arg(filetype)
  if( filetype == "csv" ){
    do_save = function(x, thing) {
      file = file.path(savepath, paste0( thing, ".", filetype ) )
      x %>% as.matrix %>% write.csv(file, quote = F, row.names = F)
    }
  } else if( filetype == "h5" ){
    stop("Sorry, hdf5 saving is not implemented yet.\n")
  } else {
    stop("Filetype must be 'csv' or 'h5' ")
  }
  looks$groups       %>% sapply(paste, collapse = " ") %>% do_save( "groups" )
  looks$vars_to_omit %>% do_save( "vars_to_omit" )
  looks$knockoffs    %>% do_save( "knockoffs" )
  for( field in names( looks$updates[[1]] ) ){
    looks$updates %>% lapply( magrittr::extract2, field ) %>% lapply(c) %>% Reduce( f = cbind ) %>% do_save( field )
  }
}

#' Deprecated. See saveLoad__loadCompactLooks
#'
#' @export
readCompactLooks = loadCompactLooks = function(...){ saveLoad__loadCompactLooks(...) }

#' Load result of saveCompactLooks
#'
#' @export
saveLoad__loadCompactLooks = function(savepath, filetype = c("csv", "h5")){
  # Parse input; set up loaders and slots
  filetype = match.arg(filetype)
  looks = list(knockoffs = NA, groups = NA, vars_to_omit = NA, updates = list(NA))
  if( filetype == "csv" ){
    do_load = function( thing ) {
      file = file.path( savepath, paste0( thing, ".", filetype ) )
      as.matrix(read.csv(file, stringsAsFactors = F))
    }
  } else if( filetype == "h5" ){
    stop("Not implemented yet.\n")
  } else {
    stop("Filetype must be 'csv' or 'h5' ")
  }
  # Get stuff
  looks$vars_to_omit = do_load( "vars_to_omit")[,1, drop = T]
  looks$knockoffs = do_load( "knockoffs" )
  # Handle irregular shape of groups and possible NULL specification of groups
  looks$groups = do_load( "groups") %>% as.character %>% strsplit(" ") %>% lapply(as.numeric)
  if(length(looks$groups)==0){
    looks$groups = NULL
  }
  # Get updates as matrices
  temp = list()
  for( field in c(
    "mean_update_left",
    "mean_update_right",
    "sqrt_cov_update",
    "random_update"
  )){
    temp[[field]] = do_load( field )
  }
  # Reshape updates into the admittedly eccentric original format
  fix_transpose = function(named_updates){
    do_transpose = grepl("right|cov", names(named_updates))
    named_updates[do_transpose] = lapply(named_updates[do_transpose], t)
    named_updates
  }
  for(k in looks$vars_to_omit){
    looks$updates[[k]] =
      lapply(temp, magrittr::extract, , k,  drop = F) %>%
      lapply(as.matrix) %>%
      fix_transpose %>%
      magrittr::set_names(names(temp))
  }
  looks
}

#' Deprecated. See saveLoad__formAllLooks
#'
#' @export
formAllLooks = function(...){ saveLoad__formAllLooks(...) }

#' Given the low-rank representations, update knockoffs to omit each variable.
#'
#' @param k variable to omit.
#' @param statistic Optional but perhaps useful for memory efficiency: instead of returning all knockoffs, compute statistics using this function, and return those instead.
#' @param updates @param knockoffs @param vars_to_omit
#' Inputs should be previously saved by \code{saveCompactLooks} and read by \code{loadCompactLooks} or from \code{generateLooks(..., output_type = 'knockoffs_compact'}.)
#' Those functions return a list with the same names as the necessary args.
#' @export
#'
saveLoad__formAllLooks = function(knockoffs, vars_to_omit, updates, statistic = NULL, X = NULL, ...){
  if(is.null(statistic)){
    return( lapply( vars_to_omit, function(k) formOneLook(knockoffs, vars_to_omit, updates, k) ) )
  } else {
    stopifnot("Pass original data to formAllLooks if you want test statistics as output.\n"=!is.null(X))
    return( lapply( vars_to_omit, function(k) statistic( X[,-k], formOneLook(knockoffs, vars_to_omit, updates, k), y = X[,k], ... ) ) )
  }
}


#' Deprecated. See saveLoad__formOneLook
#'
#' @export
formOneLook = function(...){ saveLoad__formOneLook(...) }

#' Given the low-rank representations, update knockoffs to omit one variable.
#'
#' @param k variable to omit.
#' @param updates @param knockoffs @param vars_to_omit
#' Inputs should be from \code{loadCompactLooks} or from \code{generateLooks(..., output_type = 'knockoffs_compact'}.)
#' Those functions return a list with the same names as the necessary args.
#' @export
#'
saveLoad__formOneLook = function(knockoffs, vars_to_omit, updates, k){
  one_update = getUpdateK(k, vars_to_omit, updates)
  with(one_update,
       knockoffs[,-k] +
         getMeanUpdate(one_update) +
         getRandomUpdate(one_update, n_obs = nrow(knockoffs))
  )
}
