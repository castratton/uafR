#'Custom Output from `categorate()`
#'
#'A list containing all of the output from `categorate(standard_data$Compound.Name[standard_data$Match.Factor > 78])`. Has all necessary elements for downstream functions, specifically `mzExacto()`.
#'
#'@format ## `standard_categorated`
#'A list with 4 data frames.
#'\describe{
#'  \item{Databases}{Results from Database Searches}
#'  \item{FMCS}{Results from Atomic Structure Summaries}
#'  \item{FunctionalGroups}{Results from Structural Similarity Matches with Input Library}
#'  \item{BestChemMatch}{Library Chemicals that had Strong Structural Matches (Tanimoto Similarity > 0.95)}
#'  }
"standard_categorated"
