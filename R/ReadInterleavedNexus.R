ReadInterleavedNexus <- function(nexusFile){
  dat <- scan(nexusFile, what="character", sep="\n")
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