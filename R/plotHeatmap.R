plotHeatmap <- function(SNPdataset, col=heat.colors(256)) {
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  SNPdataset <- as.matrix(SNPdataset)
  heatmatrix <- apply(SNPdataset, c(1,2), PercentMissing)
  hacked.heatmap(heatmatrix, Rowv=NA, Colv=NA, xlab="", ylab="", margins=c(5,2))

}