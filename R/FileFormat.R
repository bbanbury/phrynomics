FileFormat <- function(file){
  # returns "phy" or "nex" depending on the file type
  format <- NA
  if(any(grep("Error", try(read.table(file, skip=GetLinesToSkip(file)), silent=TRUE), ignore.case=TRUE))) {
    if(!any(grep("Error", try(read.nexus.data(file)), ignore.case=TRUE)))
      format <- "nex"
  } 
  if(!any(grep("Error", try(read.table(file, skip=GetLinesToSkip(file)), silent=TRUE), ignore.case=TRUE)))
    format <- "phy"
  return(format)
}
