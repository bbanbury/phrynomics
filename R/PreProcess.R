PreProcess <- function(SNPdataset) {
#get rid of any underscores in data (weird output from PyRad)
  if(any(which(SNPdataset == "_"))){
    underscores <- which(SNPdataset == "_", arr.ind=T)
    SNPdataset  <- SNPdataset[,-unique(underscores[,2])]
  }
  else return(SNPdataset)
}
