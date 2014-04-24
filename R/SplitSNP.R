SplitSNP <- function(SNPdataset){
#splits a single element into many (ex "ATA" to "A" "T" "A")
  loci <- dim(SNPdataset)[2]
  initialLociLengths <- nchar(SNPdataset[1,])
  splitSNP <- matrix(nrow=dim(SNPdataset)[1], ncol=sum(initialLociLengths)+loci-1)
  for(j in sequence(dim(SNPdataset)[1])) {
    splitSNP[j,] <- strsplit(paste(SNPdataset[j,], collapse=" "), "")[[1]]
  }
  rownames(splitSNP) <- rownames(SNPdataset)
  return(splitSNP)
}
