library("magrittr")
set_dimnames = function(X, dn){
  dimnames(X) = dn
  X
}

# Generate some average-joe looks
set.seed(0)
X = matrix(rnorm(1e3), ncol = 10)
looks = rlookc::generateLooks( X = X,
                               mu = 0,
                               Sigma = cor(X),
                               output_type = "knockoffs_compact",
                               method = "group",
                               groups = list(1:5, 6:10))
looks_realized = formAllLooks(knockoffs = looks$knockoffs, vars_to_omit = looks$vars_to_omit, updates = looks$updates)
dir.create("tests/test_save_looks_full/", recursive = T, showWarnings = F)
mapply(write.csv, looks_realized, paste0("tests/test_save_looks_full/", looks$vars_to_omit, ".csv"), MoreArgs = list(quote = F, row.names = F))

test_that("saveCompactLooks saves files", {
  saveCompactLooks(looks, "tests/test_save_looks")
  expect_true(file.exists( "tests/test_save_looks/knockoffs.csv" ))
  looks2 = loadCompactLooks( "tests/test_save_looks")
})

test_that("loadCompactLooks loads files correctly", {
  looks2 = loadCompactLooks( "tests/test_save_looks" )

  # Same knockoffs, vars, groups, except maybe array names
  expect_equal(looks2$knockoffs %>% set_dimnames(NULL),
               looks$knockoffs %>% set_dimnames(NULL))
  expect_equal(looks2$vars_to_omit %>% set_names(NULL),
               looks$vars_to_omit %>% set_names(NULL))
  expect_equal(looks2$groups,
               looks$groups)

  # Same updates, except maybe array names
  for(i in seq_along(looks$vars_to_omit)){
    u1 = looks$updates[[looks$vars_to_omit[[i]]]]
    u2 = looks2$updates[[looks2$vars_to_omit[[i]]]]
    expect_equal(names(u1), names(u2))
    expect_equal(u1 %>% lapply(set_dimnames, NULL),
                 u2 %>% lapply(set_dimnames, NULL))
    expect_equal(u1 %>% sapply(dim),
                 u2 %>% sapply(dim))
    expect_equal(u1 %>% lapply(c) %>% Reduce(f=c),
                 u2 %>% lapply(c) %>% Reduce(f=c))
  }
})


test_that("Loaded compact looks are realized into the correct non-compact form", {
  looks2 = loadCompactLooks( "tests/test_save_looks")
  ko = formAllLooks(looks2$knockoffs, looks2$vars_to_omit, looks2$updates)
  lapply(ko, dim) %>%  lapply(expect_equal, c(100, 9))
  mapply(
    expect_equal,
    looks_realized %>% lapply(set_dimnames, NULL),
    ko %>% lapply(set_dimnames, NULL)
  )
})
