% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardifyIt.R
\name{standardifyIt}
\alias{standardifyIt}
\title{standardifyIt}
\usage{
standardifyIt(
  data_in,
  standard_type = "Internal",
  standard_used = "Tetradecane",
  IS_ng = 190.5,
  IS_uL = 1,
  collect_time = 1,
  sample_amt = 1,
  ES_calibration = NA
)
}
\arguments{
\item{data_in}{mzExacto output}

\item{standard_type}{specifies type of standardization to perform ("Internal" or "External")}

\item{standard_used}{specifies the standard used as an internal standard (name your chemical)}

\item{IS_ng}{specifies the quantified number of molecules in the standard}

\item{IS_uL}{specifies the amount of standard added to the samples}

\item{collect_time}{specifies how long samples were collected for (e.g. hour(s), day(s), year(s))}

\item{sample_amt}{specifies how many individuals/items the samples were collected from}

\item{ES_calibration}{input matrix from which calibrations will be fit for external standardization}
}
\value{
Returns the original data standardized relative to the user-specified internal or external
inputs
}
\description{
Post-merge function to quantify compound emission rates relative to internal or external
standard(s). If using internal standard (IS), will use the contained `standardify()` function with
user-specified inputs. If using external standard (ES), an input "matrix" from which standard curves
can be derived.
}
\examples{
standardifyIt(standard_exacto, standard_type = "Internal", standard_used = "Octanal",
IS_ng = 1, IS_uL = 1, collect_time = 1)
standardifyIt(standard_exacto, standard_type = "External", ES_calibration = ExternalStandard_data)

}
