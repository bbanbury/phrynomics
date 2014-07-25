#' plot.snp
#' 
#' Plot SNP Matrix data
#' 
#' This function will create a three-panel plot, the first that shows the number of sites per locus, the second that shows the frequency distribution of sites per locus, and the third shows the amount of missing data . This doesn't work when there is a single locus, however you can plot the third panel alone by calling \code{plotMissing}.  
#' 
#' @param x an object in the class "snp"
#' @param calc Calculate missing data either by individual "sites" or "loci"
#' @param ... additional arguments passed to \code{plot}
#' @export
#' @seealso \link{ReadSNP} \link{summary.snp} \link{plotMissing}
#' @examples
#' data(fakeData)
#' myData <- ReadSNP(fakeData)
#' plot(myData)

plot.snp <- function(x, calc="loci", ...){
  dev.new(width=10, height=5)
  layout(matrix(1:3, nrow=1, byrow=TRUE), respect=TRUE)

  ## panel 1 ##
  plot(x$nsites, type="n", ylab="Number Sites", xlab="Locus", xaxt="n")
  text(1:length(x$nsites), x$nsites, labels=1:length(x$nsites), col="blue")
  segments(0, mean(x$nsites), length(x$nsites)+1, mean(x$nsites), col="red")
  text(0.1+(0.1*length(x$nsites)), mean(x$nsites)+(.02*max(x$nsites)), "mean", col="red")
  title(main="Number of Sites Per Locus")

  ## panel 2 ##
  plot(density(x$nsites), xlab="Number Sites", ylab="Frequency", main="", col="blue", ...)
  segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
  title(main="Frequency of Sites per Locus")

  ## panel 3 ##
  numSpecies <- MakePresentAbsent(x, calc=calc, returnData="presabsentdata")
  sums <- apply(numSpecies, 2, sum)
  spec <- sequence(dim(numSpecies)[1])
  y <- sapply(spec, function(z) length(which(sums >= z))) 
  plot(spec, y, pch=16, col="blue", xlab="Minimum Number of Individuals", ylab=paste("Number of", calc, "Retained"))
  lines(spec, y, col="blue")
  title(main="Missing Data")
}
























