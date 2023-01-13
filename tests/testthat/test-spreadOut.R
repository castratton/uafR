
standard_spread_t1 = suppressWarnings(spreadOut(standard_data))

test_that("outputs in list are correct sizes",{
 expect_equal(nrow(standard_spread_t1$Area), nrow(standard_data))
 expect_equal(nrow(standard_spread_t1$Compounds), nrow(standard_data))
 expect_equal(nrow(standard_spread_t1$MZ), nrow(standard_data))
 expect_equal(nrow(standard_spread_t1$MatchFactor), nrow(standard_data))
 expect_equal(nrow(standard_spread_t1$RT), nrow(standard_data))
 expect_equal(nrow(standard_spread_t1$Mass), nrow(standard_data))
 expect_equal(nrow(standard_spread_t1$rtBYmass), nrow(standard_data))
 expect_equal(length(standard_spread_t1$webInfo), nrow(standard_data))
})

test_that("items are in the correct order",{
 expect_equal(paste0(standard_spread_t1$Area[3,][!is.na(standard_spread_t1$Area[3,])]), paste0(standard_data$Component.Area[3]))
 expect_equal(paste0(standard_spread_t1$Compounds[3,][!is.na(standard_spread_t1$Compounds[3,])]), paste0(standard_data$Compound.Name[3]))
 expect_equal(paste0(standard_spread_t1$MZ[3,][!is.na(standard_spread_t1$MZ[3,])]), paste0(standard_data$Base.Peak.MZ[3]))
 expect_equal(paste0(standard_spread_t1$MatchFactor[3,][!is.na(standard_spread_t1$MatchFactor[3,])]), paste0(standard_data$Match.Factor[3]))
 expect_equal(paste0(standard_spread_t1$RT[3,][!is.na(standard_spread_t1$RT[3,])]), paste0(standard_data$Component.RT[3]))
 expect_equal(paste0(names(standard_spread_t1$webInfo[3])), paste0(standard_data$Compound.Name[3]))

 expect_equal(paste0(standard_spread_t1$Area[34,][!is.na(standard_spread_t1$Area[34,])]), paste0(standard_data$Component.Area[34]))
 expect_equal(paste0(standard_spread_t1$Compounds[34,][!is.na(standard_spread_t1$Compounds[34,])]), paste0(standard_data$Compound.Name[34]))
 expect_equal(paste0(standard_spread_t1$MZ[34,][!is.na(standard_spread_t1$MZ[34,])]), paste0(standard_data$Base.Peak.MZ[34]))
 expect_equal(paste0(standard_spread_t1$MatchFactor[34,][!is.na(standard_spread_t1$MatchFactor[34,])]), paste0(standard_data$Match.Factor[34]))
 expect_equal(paste0(standard_spread_t1$RT[34,][!is.na(standard_spread_t1$RT[34,])]), paste0(standard_data$Component.RT[34]))
 expect_equal(paste0(names(standard_spread_t1$webInfo[34])), paste0(standard_data$Compound.Name[34]))

 expect_equal(paste0(standard_spread_t1$Area[56,][!is.na(standard_spread_t1$Area[56,])]), paste0(standard_data$Component.Area[56]))
 expect_equal(paste0(standard_spread_t1$Compounds[56,][!is.na(standard_spread_t1$Compounds[56,])]), paste0(standard_data$Compound.Name[56]))
 expect_equal(paste0(standard_spread_t1$MZ[56,][!is.na(standard_spread_t1$MZ[56,])]), paste0(standard_data$Base.Peak.MZ[56]))
 expect_equal(paste0(standard_spread_t1$MatchFactor[56,][!is.na(standard_spread_t1$MatchFactor[56,])]), paste0(standard_data$Match.Factor[56]))
 expect_equal(paste0(standard_spread_t1$RT[56,][!is.na(standard_spread_t1$RT[56,])]), paste0(standard_data$Component.RT[56]))
 expect_equal(paste0(names(standard_spread_t1$webInfo[56])), paste0(standard_data$Compound.Name[56]))
})

test_that("having missing input columns is bad",{
 expect_error(spreadOut(standard_data[,-1:5]))
 expect_error(spreadOut(standard_data[,-2:6]))
 expect_error(spreadOut(standard_data[,-3:7]))
})
