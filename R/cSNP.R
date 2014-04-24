cSNP <- function(splitSNP, KeepVector=NULL, maintainLoci=TRUE){
#this function will concatenate split datasets, and reduce them either by a KeepVector and/or remove spaces indicating loci.  
  if(is.null(KeepVector))
    KeepVector <- rep(TRUE, dim(splitSNP)[2])
  catSNP <- matrix(nrow=dim(splitSNP)[1], ncol=length(which(splitSNP[1,] == " "))+1) 
  rownames(catSNP) <- rownames(splitSNP)
  for(j in sequence(dim(catSNP)[1])) {
    #totally annoying strsplit issue where it cleaves off one space at the end
    string <- strsplit(paste(splitSNP[j,][which(KeepVector=="TRUE")], collapse=""), " ")[[1]]
    if(length(string) +1 == dim(catSNP)[2]){
      string <- c(string, "")
    }
    catSNP[j,] <- string
  }
  if(any(which(catSNP[1,] == "")))
    catSNP <- catSNP[,-which(catSNP[1,] == "")]  #rewrite with lost loci
  if(!maintainLoci)
    catSNP <- apply(catSNP, 1, paste, collapse="")
  return(data.frame(catSNP, stringsAsFactors=FALSE))
}
