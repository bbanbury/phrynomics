#' Read SNP
#' 
#' Read a SNP File
#'
#' This function reads in SNP datasets into R from nexus or phylip formatted files or converts an R object dataset into one of the class "snp".  
#'
#' @aliases snp
#' @param file A file name specified by either a variable or a double-quoted string (ie, file=""). Also accepts R objects that are data frames or matrices and converts to the class "snp".  
#' @param row.names A numeric indicating which column of data the row names occur (default is 1)
#' @param preprocess If preprocess=TRUE, it will remove underscores from pyrad output ("_")
#' @export
#' @return Returns an object in the class "snp" with the following components:
#' \itemize{
#' \item \code{data} The SNP data matrix
#' \item \code{ntax} Number of total taxa present in data matrix
#' \item \code{nloci} Number of loci present in the data matrix
#' \item \code{nsnps} A vector of the number of SNPs present in each locus
#' }
#' @seealso \link{WriteSNP} 
#' @examples
#' #Read in a SNP matrix from file
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' write.table(a, file="ex.data", col.names=FALSE)
#' ReadSNP("ex.data")
#' 
#' #Read in a SNP matrix from an R object
#' data(fakeData)
#' ReadSNP(fakeData)

ReadSNP <- function(file, row.names=1, preprocess=TRUE){
  if(class(file) == "character"){
    inputFileType <- FileFormat(file)
    if(inputFileType == "phy")
      initializeTable <- read.table(file, row.names=row.names, skip=GetLinesToSkip(file), stringsAsFactors=FALSE)
    if(inputFileType == "nex") 
      initializeTable <- ConvertNexDataToPhyData(read.nexus.data(file))
    if (is.na(inputFileType))
      stop("Having a hard time reading the file")
  }
  if(class(file) == "data.frame" || class(file) == "matrix"){ 
    initializeTable <- file
  }
  colnames(initializeTable) <- paste("locus", 1:ncol(initializeTable), sep="")
  if(preprocess){
    if(any(which(initializeTable == "_"))){
      underscores <- which(initializeTable == "_", arr.ind=T)
      initializeTable <- initializeTable[,-unique(underscores[,2])]
    }
  }
  ntax <- dim(initializeTable)[1]
  nloci <- dim(initializeTable)[2]
  nsnps <- nchar(initializeTable[1,])
  snpdata <- list(data=initializeTable, ntax=ntax, nloci=nloci, nsnps=nsnps)
  class(snpdata) <- "snp"
  return(snpdata)
}