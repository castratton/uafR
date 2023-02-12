
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

## Example Mass Spectrometry Workflows

These are basic examples of how to use core functions:

``` r
library(uafR)
## example usage for Mass Spectrometry data
input_dat = read.csv("your/gcms/dataset.csv")
```
### In this example, the user knows what chemicals they are interested in:
``` r
input_spread = spreadOut(input_dat)
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")

### extract the query_chemicals from the "spread out" input:
input_exacto = mzExacto(input_spread, query_chemicals)
```
### In this example, the user just wants to keep the top hits:
``` r
query_chemicals = input_dat$Compound.Name[input_dat$Match.Factor > 80]

input_exacto = mzExacto(input_spread, query_chemicals)
```

## Example Cheminformatics Workflow
``` r
## example usage for chemical informatics:
query_chemicals = c("Linalool", "Methyl Salicylate", "Limonene", "alpha-Thujene")
GroupA = c("Guaiacol",	"Tridecane",	"Ethyl heptanoate", "Caffeine")
GroupB = c("2-Aminothiazole", "Aspirin", "Octanoic acid", "alpha-Pinene", "Toluene")
chem_library = data.frame(cbind(GroupA, GroupB))

query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")
```
## Combined Mass Spectrometry + Cheminformatics Workflow

``` r
query_chemicals = input_dat$Compound.Name[input_dat$Match.Factor > 70]
query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")

## example of using the info from categorate() to get a user-defined set of chemicals with exactoThese():
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "All")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "reactives")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "KEGG")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FEMA")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
these_chems = exactoThese(query_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Greater Than", subset_input = 125)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Less Than", subset_input = 205)
these_chems = exactoThese(query_categorated, subsetBy = "FMCS", subsetArgs = "MW", subsetArgs2 = "Between", subset_input = c(125, 200))
these_chems = exactoThese(query_categorated, subsetBy = "Library", subsetArgs = "GroupB")

input_exacto = mzExacto(input_spread, these_chems)
```
