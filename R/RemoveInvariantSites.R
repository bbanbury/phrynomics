RemoveInvariantSites <- function(SNPdataset, chatty=FALSE){
#This function will remove sites with less than two bases
#one base plus ambiguity code (that icludes base) will be dropped 
  SNPdataset <- PreProcess(SNPdataset)
  loci <- dim(SNPdataset)[2]
  initialLociLengths <- nchar(SNPdataset[1,])
  splitdata <- SplitSNP(SNPdataset)
  KeepVector <- apply(splitdata, 2, IsVariable)  
  breaks <- which(splitdata[1,] == " ")
  newSNPdataset <- cSNP(splitdata, KeepVector=KeepVector, maintainLoci=TRUE)
  newLoci <- length(which(newSNPdataset[1,] != ""))  #number of new loci 
  if(chatty)
    print(paste("removed", loci-newLoci, "of", loci, "loci"))
  return(newSNPdataset)
}
