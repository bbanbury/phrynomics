#' summary.snp
#'
#' Print Summary of SNP Matrix 
#'
#' This function prints a short summary of the SNP matrix to the console. 
#'
#' @param object an object in the class "snp"
#' @param ... additional arguments passed to \code{plot}
#' @export
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{print.snp} 
#' @examples
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' write.table(a, file="ex.data", col.names=FALSE)
#' b <- ReadSNP("ex.data")
#' summary(b)
#' str(b)

summary.snp <- function (object, ...){
  cat("\nSNP dataset:", deparse(substitute(object)), "\n")
  cat(paste("\nSNP dataset contains",  object$ntax, "taxa,", object$nloci, "loci, and", sum(object$nsites), "sites\n\n"))
  cat("\tNumber of taxa:", object$ntax, "\n")
  cat("\tNumber of loci:", object$nloci, "\n")
  cat("\tNumber of sites:", sum(object$nsites), "\n")
  cat("\n\tAverage number of sites per locus:", mean(object$nsites), "\n")
  cat("\tMinimum number of sites per locus:", min(object $nsites), "\n")
  cat("\tMaximum number of sites per locus:", max(object $nsites), "\n")
}