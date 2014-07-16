#' Fake Data
#' 
#' A fake dataset that includes three loci with 10 SNPs with various chatracteristics.  We use this dataset in our examples to illustrate how some of the phrynomics functions work to subset data matrices.  Locus 1 has four SNPs, Loci 2 and 3 have three SNPs. These SNPs were created with the following in mind:
#' \tabular{lll}{
#' Locus 1: \tab "A" and missing \tab Should fail binary and variable check. \cr
#' Locus 2: \tab "A", "T", and missing \tab Should be good. \cr
#' Locus 3: \tab "A", "T", "G", and missing data \tab Should fail binary check. \cr
#' Locus 4: \tab "A", "T", "R", and missing data. \tab Because "R" is an ambiguity code for "A" and "T", this is both binary and variable. \cr
#' Locus 5: \tab "A" and "R" \tab Is binary but not variable, as the ambiguity code includes "A".\cr
#' Locus 6: \tab "A", "D", and "R" \tab Is not Binary, because "D" is ambiguous for three alleles. Is not variable, as both ambiguity codes include "A". \cr
#' Locus 7: \tab "A" and missing \tab Should fail binary and variable check. \cr
#' Locus 8: \tab "A" and missing \tab Should fail binary and variable check. \cr
#' Locus 9: \tab "A", "G", and missing \tab Is variable and binary. \cr
#' Locus 10: \tab "A", "G", "R" \tab Is Variable and binary. \cr
#' }
#' docType data
#' @keywords datasets
#' @name fakeData
NULL