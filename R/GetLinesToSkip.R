GetLinesToSkip <- function(file){
  a <- length(suppressWarnings(system(paste("grep 'ID:' ", file, sep=""), intern=TRUE)))
  b <- length(suppressWarnings(system(paste("grep ^[#+] ", file, sep=""), intern=TRUE)))
  return(a+b)
}