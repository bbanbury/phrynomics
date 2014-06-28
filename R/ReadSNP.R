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