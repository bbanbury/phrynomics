#' File Format 
#' 
#' This is a function to autodetect phylip or nexus style file formatting.
#' 
#' @param file A file name specified by either a variable or a double-quoted string (ie, file="")
#' @export
#' @return Returns either "phy" or "nex"
#' @seealso \link{ReadSNP} 
#' @examples
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' WriteSNP(a, "example.phy")
#' FileFormat("example.phy")
#'
#' WriteSNP(a, "example.nex", format="nexus") 
#' FileFormat("example.nex")

FileFormat <- function(file){
  format <- NA
  if(any(grep("Error", try(read.table(file, skip=GetLinesToSkip(file)), silent=TRUE), ignore.case=TRUE))) {
    if(!any(grep("Error", try(read.nexus.data(file)), ignore.case=TRUE)))
      format <- "nex"
  } 
  if(!any(grep("Error", try(read.table(file, skip=GetLinesToSkip(file)), silent=TRUE), ignore.case=TRUE)))
    format <- "phy"
  return(format)
}
