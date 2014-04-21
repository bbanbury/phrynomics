#Get taxa data overlap

dataset <- read.table("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/c70p3noAmbigs.snps", row.names=1, colClasses="character", skip=1)
dataset <- read.table("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/c5p3noAmbigs.snps", row.names=1, colClasses="character", skip=1)

source("~/phrynomics/trunk/phrynomicsFunctions.R")


MakePresentAbsent <- function(dataset){
#function removes loci breaks, and changes all data to 1, all missing data to 0
  split <- SplitSNP(dataset)
  split <- split[,-which(split[1,] == " ")]
  split[c(which(split=="-"), which(split=="N"))] <- 0
  split[-c(which(split=="-"), which(split=="N"), which(split=="0"))] <- 1
  return(matrix(as.numeric(unlist(split)),nrow=nrow(split)))
}


DataOverlap <- function (dataset){
#This function will calculate the amount of shared SNPs
  presAbsentDataset <- MakePresentAbsent(dataset)
  rownames(presAbsentDataset) <- rownames(dataset)
  totals <- apply(presAbsentDataset, 2, sum)
  probs <- totals/73
  for(i in sequence(dim(presAbsentDataset)[2])){  #probably could clean this up by using replace()
    presAbsentDataset[which(presAbsentDataset[,i] == 1),i] <- probs[i]
  }  
  speciesTotals <- apply(presAbsentDataset, 1, mean)
  return(list(presAbsentDataset, speciesTotals))
}
DataOverlap(dataset)[[2]]



#calculate branch length error from root to tips
#use BL.AllTrees
BL.AllTrees <- BL.AllTrees[[1]]

GetAncestors <- function(treeEdgeMatrix, tip) {
#function to tree traverse back to the root to get node numbers
  Root <- min(treeEdgeMatrix[,1])
  is.done <- FALSE
  desc <- tip
  desc.vector <- tip
  while(!is.done){
    a <- which(treeEdgeMatrix[, 2] == desc)
    b <- treeEdgeMatrix[a, 1]
    desc.vector <- c(desc.vector, b)
    if(b == Root)
      is.done <- TRUE
    else
      desc <- b
  }
  return(desc.vector)
}


CalculateTotalTipBLError <- function(BL.AllTrees) {
#this function should take BL.AllTrees matrix and return root to tip totals BL.DIFF
  tips <- BL.AllTrees[which(BL.AllTrees$class == "tip"), 2]
  tipPathDifferences <- NULL
  for(i in tips){
    ancs <- GetAncestors(BL.AllTrees[,1:2], tips[i])
    pathDifferences <- BL.AllTrees[which(BL.AllTrees[,2] %in% ancs), 12]  #column 12 is relativeBLdiff (which can be pos or neg)
    totalPathDifference <- sum(abs(pathDifferences))  # sum absolute value path differences # relative number
    tipPathDifferences <- c(tipPathDifferences, totalPathDifference)
  }
  names(tipPathDifferences) <- paste("tip", tips, sep="")
  return(tipPathDifferences)
}

CalculateTotalTipBLError(BL.AllTrees)

GetJustTipBLError <- function(BL.AllTrees){
#this function will return a vector of BL differences for just tips
  tips <- BL.AllTrees[which(BL.AllTrees$class == "tip"), 2]
  tipDifferences <- rep(0, length(tips))
  for(i in tips){
    tipDifferences[i] <- BL.AllTrees[which(BL.AllTrees[,2] == tips[i]), 12]  #column 12 is relativeBLdiff (which can be pos or neg)
    names(tipDifferences)[i] <- paste("tip", tips[i], sep="")
  }
  return(tipDifferences)
}

GetJustTipBLError(BL.AllTrees)






















