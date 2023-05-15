#' plot.snp
#' 
#' Plot SNP Matrix data
#' 
#' This function will create a plot of the amount of missing data. It is calculated as the number of retained loci or sites given the minimum taxa that must be included.  For example, if you want a dataset that only includes sites with all your individuals (ie NO missing data), you can see how many sites will be included.
#' 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param calc Calculate missing data either by individual "sites" or "loci"
#' @param ... additional arguments passed to \code{plot}
#' @export
#' @seealso \link{ReadSNP} \link{summary.snp} \link{plot.snp}
#' @examples
#' data(fakeData)
#' myData <- ReadSNP(fakeData)
#' plotMissing(myData, calc="loci")
#' plotMissing(myData, calc="sites")

plotMissing <- function(SNPdataset, calc="loci", ...){
  if(inherits(SNPdataset, "snp"))
    SNPdataset <- SNPdataset$data
  numSpecies <- MakePresentAbsent(SNPdataset, calc=calc)
  sums <- apply(numSpecies, 2, sum)
  x <- sequence(dim(numSpecies)[1])
  y <- sapply(x, function(x) length(which(sums >= x))) 
  plot(x, y, pch=16, col="blue", xlab="Minimum Number of Individuals", ylab=paste("Number of", calc, "Retained"), ...)
  lines(x, y, col="blue")
  title(main="Missing Data")
}