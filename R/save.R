#' Save all the leave-one-out knockoffs to disk, e.g. for reading into Python.
#'
#' Save them in their low-rank representation to avoid wasting space.
#'
saveCompactLooks = function(looks, savepath, filetype = c("csv", "h5")){
  dir.create(savepath, recursive = T, showWarnings = F)
  filetype = match.arg(filetype)
  if( filetype == "csv" ){
    do_save = function(x, thing) {
      file = file.path(savepath, paste0( thing, ".", filetype ) )
      x %>% as.matrix %>% write.csv(file, quote = F, row.names = F)
    }
  } else if( filetype == "h5" ){
    stop("Not implemented yet.\n")
  } else {
    stop("Filetype must be 'csv' or 'h5' ")
  }
  looks$groups       %>% sapply(paste, collapse = " ") %>% do_save( "groups" )
  looks$vars_to_omit %>% do_save( "vars_to_omit" )
  looks$knockoffs    %>% do_save( "knockoffs" )
  for( field in names( looks$updates[[1]] ) ){
    looks$updates %>% lapply( extract2, field ) %>% lapply(c) %>% Reduce( f = cbind ) %>% do_save( field )
  }
}

#' Load result of saveCompactLooks
#'
loadCompactLooks = function(savepath, filetype = c("csv", "h5")){
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
    "mean_update1_left",
    "mean_update1_right",
    "mean_update2_left",
    "mean_update2_right",
    "sqrt_cov_update",
    "random"
  )){
    temp[[field]] = do_load( field )
  }
  # Reshape updates into the admittedly eccentric original format
  fix_transpose = function(named_updates){
    do_transpose = grepl("right|cov", names(named_updates))
    named_updates[do_transpose] %<>% lapply(t)
    named_updates
  }
  for(k in looks$vars_to_omit){
    looks$updates[[k]] = lapply(temp, extract, , k,  drop = F) %>% lapply(as.matrix) %>% fix_transpose %>% set_names(names(temp))
  }
  looks
}

#' Alias for loadCompactLooks
#'
readCompactLooks = loadCompactLooks
