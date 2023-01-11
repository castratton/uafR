#'Standard GC/MS Output
#'
#'A subset of data from a study running 3 samples containing
#'known quantities of standard chemicals -- (1) ethyl hexanoate,
#'(2) methyl salicylate, (3) octanal, and (4) undecane.
#'
#'@format ## `standard_data`
#'A data frame with 63 rows and 8 columns:
#'\describe{
#'  \item{Component.RT}{Retention Time}
#'  \item{Base.Peak.MZ}{Captured M/Z Value}
#'  \item{Base.Peak.Area}{Area/Quantity of Chemical}
#'  \item{Compound.Name}{Name of Tentative Chemical}
#'  \item{Match.Factor}{Accuracy of Tentative Match}
#'  \item{Sample.Name}{Sample name}
#'  \item{File.Name}{Redundant - use sample name}
#'  }
"standard_data"
