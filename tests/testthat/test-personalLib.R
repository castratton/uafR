
test_that("it works!", {
 expect_no_message(personalLib(library_data, "wide", "vectors"))
 expect_no_message(personalLib(library_data, "wide", "list"))
 expect_no_message(personalLib(library_data_long, "long", "vectors"))
 expect_no_message(personalLib(library_data_long, "long", "list"))
})

test_that("long doesn't work with the wrong input", {
 expect_error(personalLib(library_data, "long", "vectors"))
 expect_error(personalLib(library_data, "long", "list"))
})
