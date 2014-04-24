ConvertMissingData <- function(SNPdataset, oldmissing="N", newmissing="?"){
#This function will gsub oldMissing to newMissing
  for(i in sequence(dim(SNPdataset)[1])){
    SNPdataset[i,] <- gsub(oldmissing, newmissing, SNPdataset[i,])
  }
  return(SNPdataset)
}
