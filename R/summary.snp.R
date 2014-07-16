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
  cat(paste("\nSNP dataset contains",  x$ntax, "taxa,", x$nloci, "loci, and", sum(x$nsnps), "SNPs\n\n"))
  cat("\tNumber of taxa:", x$ntax, "\n")
  cat("\tNumber of loci:", x$nloci, "\n")
  cat("\tNumber of snps:", sum(x$nsnps), "\n")
  cat("\n\tAverage number of SNPs per locus:", mean(x$nsnps), "\n")
  cat("\tMinimum number of SNPs per locus:", min(x$nsnps), "\n")
  cat("\tMaximum number of SNPs per locus:", max(x$nsnps), "\n")
}