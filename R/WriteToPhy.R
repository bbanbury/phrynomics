WriteToPhy <- function(data, fileName="file.phy") {
  file <- paste(dim(data)[1], nchar(paste(data[1,], collapse="")))
  write(file, file=fileName)
  write.table(data, file=fileName, append=TRUE, quote=FALSE, col.names=FALSE)  
}