#' Read Interleaved Nexus File
#' 
#' Read Interleaved Nexus File
#'
#' This function reads in SNP datasets that are interleaved and converts them to data frame. These can then be written as standard nexus or phylip formatted files.   
#'
#' @aliases ReadInterleavedNexus
#' @param file A file name specified by either a variable or a double-quoted string (ie, file=""). Also accepts R objects that are data frames or matrices and converts to the class "snp".  
#' @export
#' @return Returns a data frame with rownames as species and column(s) as loci
#' @seealso \link{WriteSNP} 

ReadInterleavedNexus <- function(file){
#make sure this works with multiple loci?
  dat <- scan(file, what="character", sep="\n")
  if(length(grep("\t;", dat) > 1))
    dat <- dat[-grep("\t;", dat)]
  whereDataStarts <- grep("matrix", dat, ignore.case=TRUE)+1
  whereDataEnds <- grep("end", dat)[which(grep("end", dat) > whereDataStarts)][1] -1
  dat <- dat[whereDataStarts:whereDataEnds]
  datlist <- strsplit(dat, split="[ +]", perl=TRUE)
  dattable <- matrix(nrow=length(datlist), ncol=2)
  for(i in sequence(length(datlist))){
    dattable[i,] <- c(datlist[[i]][which(datlist[[i]] != "")])
  }
  newdattable <- matrix(nrow=length(unique(dattable[,1])), ncol=1)
  rownames(newdattable) <- unique(dattable[,1])
  for(i in sequence(dim(newdattable)[1])){
    newdattable[i,1] <- paste(dattable[which(dattable[,1] == rownames(newdattable)),2], collapse="")
  }
  return(as.data.frame(newdattable))
}