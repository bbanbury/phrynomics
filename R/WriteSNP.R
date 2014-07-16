#' Write SNP
#' 
#' Write a SNP File
#'
#' This function writes SNP datasets into nexus or phylip formatted files. 
#'
#' @param data A data matrix of the class "matrix", "data.frame", or "snp"
#' @param file A file name specified by either a variable or a double-quoted string. If not specified (ie, file=""), then snp data is printed to the console. 
#' @param format Format for file, either "nexus" or "phylip"
#' @param missing A character denoting missing data, usually "N", "?", or "-"
#' @export
#' @seealso \link{ReadSNP} 
#' @examples
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' write.table(a, file="ex.data", col.names=FALSE)
#' b <- ReadSNP("ex.data")
#' WriteSNP(b, file="refakeData.phy")
#'
#' data(fakeData)
#' WriteSNP(fakeData, file="refakeData.nex", format="nexus")

WriteSNP <- function(data, file="", format="phylip", missing="N"){
  if(class(data) == "snp")
    data <- data$data
  if(format == "phylip"){
    write(c(dim(data)[1], sum(nchar(data[1,]))), file=file)
    write.table(data, file=file, quote=FALSE, append=TRUE, col.names=FALSE)
  }
  if(format == "nexus"){
    if(class(data) == "data.frame" || class(data) == "matrix")
      if(dim(data)[2] > 1)
        data <- data.frame(apply(data, 1, paste, collapse=""))
    nchars <- min(apply(data, 1, nchar))
    write(paste("#NEXUS"), file)
    write(paste("[Written ", Sys.Date(), " via phrynomics]", sep=""), file, append=TRUE)
    write(paste("BEGIN Data;"), file, append=TRUE)
    write(paste("   DIMENSIONS NTAX=", length(data), " NCHAR=", nchars, ";", sep=""), file, append=TRUE)
    write(paste("   FORMAT DATATYPE=Standard INTERLEAVE=no missing=", missing, ";", sep=""), file, append=TRUE)
    write(paste("Matrix"), file, append=TRUE)
    write.table(data, file, append=TRUE, quote=FALSE, col.names=FALSE)  
    write(paste(""), file, append=TRUE)
    write(paste(";"), file, append=TRUE)
    write(paste("END;"), file, append=TRUE)
    write(paste(""), file, append=TRUE)
  }
}