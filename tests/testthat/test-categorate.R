
query_chemicals = c("Linalool", "alpha-Pinene", "Aspirin", "Caffeine", "Limonene")
query_categorated = suppressWarnings(categorate(query_chemicals, library_data, input_format = "wide"))

test_that("output has the correct number of objects",{
  expect_equal(length(query_categorated), 4)
})

test_that("only items in output are for query chemicals",{
 expect_equal(unique(query_categorated$Databases$Chemical), query_chemicals)
 expect_equal(length(unique(query_categorated$FMCS$Chemical)), length(query_chemicals))
 expect_equal(length(unique(query_categorated$FunctionalGroups$Chemical)), length(query_chemicals))
 expect_equal(length(unique(query_categorated$BestChemMatch$Chemical)), length(query_chemicals))
})
