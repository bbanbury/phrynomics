GetBaseFrequencies <- function(SNP){
#note these may not add up to 1 if there is missing data or heteros
  bases <- ReturnUniqueBases(SNP)
  return(sapply(bases, function(x) length(which(SNP == x))/length(SNP)))
}
