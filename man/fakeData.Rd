% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fakeData.R
\name{fakeData}
\alias{fakeData}
\title{Fake Data}
\description{
A fake dataset that includes three loci with 10 sites with various chatracteristics.  We use this dataset in our examples to illustrate how some of the phrynomics functions work to subset data matrices.  Locus 1 has four sites, Loci 2 and 3 have three sites. These loci were created with the following in mind:
\tabular{lll}{
site 1: \tab "A" and missing \tab Should fail binary and variable check. \cr
site 2: \tab "A", "T", and missing \tab Should be good. \cr
site 3: \tab "A", "T", "G", and missing data \tab Should fail binary check. \cr
site 4: \tab "A", "T", "R", and missing data. \tab Because "R" is an ambiguity code for "A" and "T", this is both binary and variable. \cr
site 5: \tab "A" and "R" \tab Is binary but not variable, as the ambiguity code includes "A".\cr
site 6: \tab "A", "D", and "R" \tab Is not Binary, because "D" is ambiguous for three alleles. Is not variable, as both ambiguity codes include "A". \cr
site 7: \tab "A" and missing \tab Should fail binary and variable check. \cr
site 8: \tab "A" and missing \tab Should fail binary and variable check. \cr
site 9: \tab "A", "G", and missing \tab Is variable and binary. \cr
site 10: \tab "A", "G", "R", and "V" \tab Is Variable and binary. \cr
}
docType data
}
\keyword{datasets}
