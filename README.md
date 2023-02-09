
# uafR - A new standard for mass spectrometry data processing

<!-- badges: start -->
<!-- badges: end -->

## Objective

An R package that automates GC/LC-MS processing.

## Installation

When made public, you can install the development version of uafR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("castratton/uafR")
```

## Example GC-MS + Cheminformatics Workflows

These are basic examples of how to use core functions:

``` r
library(uafR)
## example usage for Mass Spectrometry data
input_dat = read.csv("your/gcms/dataset.csv")
input_spread = spreadOut(input_dat)

### in this example, the user knows what chemicals they are interested in:
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")

### extract these chemicals from the "spread out" input:
input_exacto = mzExacto(input_spread, query_chemicals)

## example usage for chemical informatics:
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")
GroupA = c("Guaiacol",	"Tridecane",	"Ethyl heptanoate", "Caffeine")
GroupB = c("2-Aminothiazole", "Aspirin", "Octanoic acid", "alpha-Pinene", "Toluene")
chem_library = data.frame(cbind(GroupA, GroupB))

query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")

## example of using the info from categorate() to get a user-defined set of chemicals:
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "All")
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "reactives")
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "KEGG")
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FEMA")
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
exactoThese(query_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Greater Than", subset_input = 125)
exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Less Than", subset_input = 205)
exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Between", subset_input = c(125, 200))
exactoThese(query_categorated, subsetBy = "Library", subsetArgs = "GroupB")
```
