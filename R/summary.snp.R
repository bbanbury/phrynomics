#' summary.snp
#'
#' Print Summary of SNP Matrix 
#'
#' This function prints a short summary of the SNP matrix to the console. 
#'
#' @param x an object in the class "snp"
#' @export
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{print.snp} 
#' @examples
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' write.table(a, file="ex.data", col.names=FALSE)
#' b <- ReadSNP("ex.data")
#' summary(b)
#' str(b)

summary.snp <- function (x, ...){
  cat("\nSNP dataset:", deparse(substitute(x)), "\n")
  cat(paste("\nSNP dataset contains",  x$ntax, "taxa,", x$nloci, "loci, and", sum(x$nsites), "SNPs\n\n"))
  cat("\tNumber of taxa:", x$ntax, "\n")
  cat("\tNumber of loci:", x$nloci, "\n")
  cat("\tNumber of sites:", sum(x$nsites), "\n")
  cat("\n\tAverage number of sites per locus:", mean(x$nsites), "\n")
  cat("\tMinimum number of sites per locus:", min(x$nsites), "\n")
  cat("\tMaximum number of sites per locus:", max(x$nsites), "\n")
}