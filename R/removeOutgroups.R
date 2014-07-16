removeOutgroups <- function(SNPdataset, taxon="PH"){
#probably only relevant to our study
  return(SNPdataset[grep(taxon, rownames(SNPdataset)),])
}
