#' Split SNP
#' 
#' SplitSNP breaks character strings into single character elements. For example, the string "ATA" will return "A" "T" "A". Breaks in loci are maintained.
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @export
#' @return Returns a matrix with a single SNP in each colomn. Empty columns represent breaks in loci. 
#' @seealso \link{ReadSNP} \link{cSNP}
#' @examples
#' data(fakeData)
#' fakeData <- ReadSNP(fakeData)
#' SplitSNP(fakeData)
#' SplitSNP(fakeData)$data

SplitSNP <- function(SNPdataset){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  loci <- dim(SNPdataset)[2]
  initialLociLengths <- nchar(SNPdataset[1,])
  splitSNP <- data.frame(matrix(nrow=dim(SNPdataset)[1], ncol=sum(initialLociLengths)+loci-1))
  for(j in sequence(dim(SNPdataset)[1])) {
    splitSNP[j,] <- strsplit(paste(SNPdataset[j,], collapse=" "), "")[[1]]
  }
  rownames(splitSNP) <- rownames(SNPdataset)
  if(snpclass == "snp")
    return(ReadSNP(splitSNP))
  else
    return(splitSNP)
}
