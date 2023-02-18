#'Custom Output from `mzExacto()`
#'
#'A data set containing all of the output from `mzExacto(standard_spread, query_chemicals)`
#'`mzExacto()`
#'
#'@format ## `standard_exacto`
#'A data frame with 4 rows and 7 columns.
#'\describe{
#'  \item{Compound}{Column with Chemical Names}
#'  \item{Mass}{Column with Published Chemical Masses}
#'  \item{RT}{Column with Best Identified Retention Times}
#'  \item{Best Match}{Column with the Top Match Factors}
#'  \item{Std_soln_00.D}{Column with Component Areas for First Sample}
#'  \item{Std_soln_07.D}{Column with Component Areas for Second Sample}
#'  \item{Std_soln_00a.D}{Column with Component Areas for Third Sample}
#'  }
"standard_exacto"
