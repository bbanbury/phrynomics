#' Get Lines to Skip 
#' 
#' This function will auto calculate the number of lines to skip before reading in either a phylip or nexus formatted file. It skips phylip-format summary lines and ones that are commented out or have ID tags (for example: ##, [ID: XX]).
#' 
#' @param file A file name specified by either a variable or a double-quoted string (ie, file="")
#' @export
#' @return Returns a number of lines to skip
#' @seealso \link{ReadSNP} 
#' @examples
#' a <- matrix(data=c("AANNTGG", "AATTTGC", "TAAATGC"), dimnames=list(c("A", "B", "C"), "locus"))
#' WriteSNP(a, "example.phy")
#' GetLinesToSkip("example.phy")

GetLinesToSkip <- function(file){
  a <- length(suppressWarnings(system(paste("grep 'ID:' ", file, sep=""), intern=TRUE)))
  b <- length(suppressWarnings(system(paste("grep ^[#+] ", file, sep=""), intern=TRUE)))
  c <- 0
  if(length(strsplit(system(paste("head -n 1 ", file, sep=""), intern=TRUE), split=" ")[[1]]) == 2){
    if(!system(paste("head -n 2 ", file, sep=""), intern=TRUE)[1] == system(paste("head -n 2 ", file, sep=""), intern=TRUE)[2])
      c <- 1
  }
  return(a+b+c)
}