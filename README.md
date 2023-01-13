
# uafR - A new standard for mass spectrometry data processing

<!-- badges: start -->
<!-- badges: end -->

## Objective

An R package that automates GC/LC-MS processing.

## Installation

You can install the development version of uafR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("castratton/uafR")
```

## Example

This is a basic example which shows you how to solve a common problem:

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
GroupA = c("Guaiacol",	"Tridecane",	"Ethyl heptanoate")
GroupB = c("2-Aminothiazole", "Aspirin", "Octanoic acid")
chem_library = data.frame(cbind(GroupA, GroupB))

query_categorated = categorate(query_chemicals, chem_library, input_format = "wide")
```

![](inst/images/TLI_USDA.png)
