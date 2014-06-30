#' Read SNP
#' 
#' Read a SNP File
#'
#' This function reads in SNP datasets into R from nexus or phylip formatted files. 
#'
#' @alias snp
#' @param file A file name specified by either a variable or a double-quoted string (ie, file="")
#' @param row.names A numeric indicating which column of data the row names occur (default is 1)
#' @export
#' @return Returns an object in the class "snp" with the following components:
#' @item data The SNP data matrix
#' @item ntax Number of total taxa present in data matrix
#' @item nloci Number of loci present in the data matrix
#' @item nsnps A vector of the number of SNPs present in each locus
#' @seealso \link{WriteSNP} 
#' @author B.Banbury
#' @examples
#' 
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' write.table(a, file="ex.data", col.names=FALSE)
#' ReadSNP("ex.data")

ReadSNP <- function(file, row.names=1){
  inputFileType <- FileFormat(file)
  if(inputFileType == "phy")
    initializeTable <- read.table(file, row.names=row.names, skip=GetLinesToSkip(file), stringsAsFactors=FALSE)
  if(inputFileType == "nex") 
    initializeTable <- convertNexDataToPhyData(read.nexus.data(file))
  if (is.na(inputFileType))
    stop("Having a hard time reading the file")
  colnames(initializeTable) <- paste("locus", 1:ncol(initializeTable), sep="")
  ntax <- dim(initializeTable)[1]
  nloci <- dim(initializeTable)[2]
  nsnps <- nchar(initializeTable[1,])
  snpdata <- list(data=initializeTable, ntax=ntax, nloci=nloci, nsnps=nsnps)
  class(snpdata) <- "snp"
  return(snpdata)
}