#' print.snp
#'
#' Compact Display of SNP Matrix 
#'
#' This function prints a short summary of the SNP matrix to the console.  
#'
#' @param x an object in the class "snp"
#' @export
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{summary.snp} 
#' @examples
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' write.table(a, file="ex.data", col.names=FALSE)
#' b <- ReadSNP("ex.data")
#' print(b)
#' str(b)

print.snp <- function (x, ...){
  cat(paste("\nSNP dataset contains",  x$ntax, "taxa,", x$nloci, "loci, and", sum(x$nsites), "sites\n\n"))
  cat("Taxa include:\n")
  cat(paste("\t", paste(rownames(x$data)[1:min(6, x$ntax)], collapse = ", "), ", ...\n", sep = ""))
}
