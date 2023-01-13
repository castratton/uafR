

test_that("output is always correct size", {
  search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))

  search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))

  search_chems = c("ethyl hexanoate", "undecane")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))

  search_chems = c("undecane")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))
})

test_that("duplicates do not matter", {
 search_chems = c("ethyl hexanoate", "ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))

 search_chems = c("ethyl hexanoate", "ethyl hexanoate", "methyl salicylate", "methyl salicylate", "octanal", "undecane")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane", "methyl salicylate", "octanal")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane",
                  "ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))
})

test_that("missing query chemicals are bad",{
 search_chems = c("", "ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
 expect_error(mzExacto(standard_spread, search_chems))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "", "octanal", "undecane")
 expect_error(mzExacto(standard_spread, search_chems))

 search_chems = c("", "", "methyl salicylate", "octanal", "undecane")
 expect_error(mzExacto(standard_spread, search_chems))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane", "", "")
 expect_error(mzExacto(standard_spread, search_chems))
})
