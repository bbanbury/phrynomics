#' Plot SNP Missing Data Heatmap
#' 
#' This function will create a heatmap plot of the level of missing data. Colors represent degrees of missing data, completely missing is opaque color.   
#' 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param calc Calculate missing data either by individual "sites" or "loci". Calculating by "loci" will give a percentage of missing data for each cell.  
#' @param col Optional color parameter that represents data (missing data will be white)  
#' @param colorMissing If set to TRUE, then missing data will be colored and white will be complete.   
#' @param xaxisTicks Number of tick marks for the x-axis   
#' @param grid Include background grid   
#' @param ... additional arguments passed to \code{plot}
#' @export
#' @seealso \link{ReadSNP} \link{CalculateMissingData} \link{PercentMissing}
#' @examples
#' data(fakeData)
#' plotHeatmap(fakeData, calc="loci")
#' plotHeatmap(fakeData, calc="sites", "goldenrod")
#' 
#' exdata <- ReadSNP(system.file("extdata", "SNPdata.snps", package="phrynomics"))
#' plotHeatmap(exdata, calc="loci")
#' plotHeatmap(exdata, colorMissing=TRUE, grid=FALSE, col="red", xaxisTicks=0, cex.axis=0.5)

plotHeatmap <- function(SNPdataset, calc="loci", col="black", colorMissing=FALSE, grid=TRUE, xaxisTicks=10, ...) {
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  if(calc == "sites"){
    SNPdataset <- SplitSNP(SNPdataset)
    if(any(SNPdataset == " "))
      SNPdataset <- SNPdataset[,-which(SNPdataset[1,] == " ")]
  }
  SNPdataset <- as.matrix(SNPdataset)
  heatmatrix <- apply(SNPdataset, c(1,2), PercentMissing)
  heatmatrix <- heatmatrix[rev(rownames(heatmatrix)),]  #reverse order of mat for plotting bottom-up
  a <- c("complete", "missing")
  if(!colorMissing){
    a <- rev(a)
    heatmatrix <- abs(heatmatrix - 1)  #so missing data gets plotted as white
  }
  ColTrans <- col2rgb(col, TRUE)/255
  alpha <- seq(0, 1, length=50)
  ColTrans <- ColTrans[, rep(1, length(alpha)), drop=FALSE]
  newcols <- rgb(ColTrans[1,], ColTrans[2,], ColTrans[3,], alpha)
  ColorLevels <- seq(0, 1, length=length(ColTrans))
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  par(mar = c(5,5,2.5,1), font=2)
  image(1:dim(heatmatrix)[2], 1:dim(heatmatrix)[1], t(heatmatrix), col=newcols, axes=FALSE, xlab=calc, ylab="Individuals", zlim=c(0,1), ...)
  box()
  if(grid)
    grid()
  if(!is.null(colnames(heatmatrix)))
    xax <- colnames(heatmatrix)
  else
    xax <- sequence(dim(heatmatrix)[2])
  xaxpos <- sequence(length(xax))
  if(length(xax) > xaxisTicks){
    xaxpos <- floor(seq.int(from=1, to=length(xax), length.out=xaxisTicks))
    xax <- xax[xaxpos]
  }
  yax <- rownames(heatmatrix)
  axis(side=1, at=xaxpos, labels=xax, cex.axis=par()$cex*0.75, ...)
  axis(side=2, at=seq(1,length(yax),1), labels=yax, las=1, cex.axis=par()$cex*0.75, ...)
  if(calc == "loci"){
    image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1), col=newcols, xlab="", ylab="", xaxt="n", yaxt="n", las=1, ...)
    ylabs <- c(a[1], seq(0.1,.9, 0.1), a[2])
    axis(side=2, at=seq(0,1, 0.1), labels=ylabs, las=1, cex.axis=par()$cex*0.75, ...)
  }
  if(calc == "sites"){
    ColorLevels <- c(ColorLevels[1], ColorLevels[length(ColorLevels)])
    image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1), col=newcols, xlab="", ylab="", xaxt="n", yaxt="n", las=1, ...)
    axis(side=2, at=c(0,1), labels=a, las=1, cex.axis=par()$cex*0.75, ...)
  }
}