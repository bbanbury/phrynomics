ReturnNucs <- function(NucCode, forSNAPP=FALSE) {
#Function for returning possible bases
#SNAPP missing data needs to be handled differently, so for now just return "-"
  possibilities <- NULL
  if(NucCode == "A" || NucCode == "G" || NucCode == "C" || NucCode == "T" || NucCode == "U")
    possibilities <- NucCode
  if(NucCode == "N" || NucCode == "-" || NucCode == "?") {
    if(forSNAPP)  possibilities <- "-"
    else  possibilities <- c("A", "G", "C", "T", "U")
  }
  if(NucCode == "R")  possibilities <- c("A", "G")
  if(NucCode == "Y")  possibilities <- c("C", "T")
  if(NucCode == "W")  possibilities <- c("A", "T")
  if(NucCode == "S")  possibilities <- c("G", "C")
  if(NucCode == "M")  possibilities <- c("A", "C")
  if(NucCode == "K")  possibilities <- c("G", "T")
  if(NucCode == "B")  possibilities <- c("G", "C", "T")
  if(NucCode == "H")  possibilities <- c("A", "C", "T")
  if(NucCode == "D")  possibilities <- c("A", "G", "T")
  if(NucCode == "V")  possibilities <- c("A", "G", "C")
  return(possibilities)
}
