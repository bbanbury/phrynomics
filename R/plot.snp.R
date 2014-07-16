#' plot.snp
#' 
#' Plot SNP Matrix data
#' 
#' This function will create a two-panel plot, the first that shows the number of SNP sites per locus and the second that shows the frequency distribution of SNP's per locus. This doesn't work when there is a single locus. 
#' 
#' @param x an object in the class "snp"
#' @param ... additional arguments passed to \code{plot}
#' @export
#' @seealso \link{ReadSNP} \link{summary.snp} 
#' @examples
#' data(fakeData)
#' snps <- ReadSNP(fakeData)
#' plot(snps)

plot.snp <- function(x, ...){
  dev.new(width=10, height=5)
  layout(matrix(1:2, nrow=1, byrow=TRUE), respect=TRUE)
  plot(x$nsnps, type="n", ylab="Number SNPs", xlab="Locus", xaxt="n")
  text(1:length(x$nsnps), x$nsnps, labels=1:length(x$nsnps), col="blue")
  segments(0, mean(x$nsnps), length(x$nsnps)+1, mean(x$nsnps), col="red")
  text(0.1+(0.1*length(x$nsnps)), mean(x$nsnps)+(.02*max(x$nsnps)), "mean", col="red")
  title(main="Number of SNPs Per Locus")
  plot(density(x$nsnps), xlab="Number SNPs", ylab="Frequency", main="", col="blue", ...)
  segments(mean(x$nsnps), -0.1, mean(x$nsnps), max(1, max(density(x$nsnps)$y)), col="red")
  title(main="Frequency of SNPs per Locus")
}
