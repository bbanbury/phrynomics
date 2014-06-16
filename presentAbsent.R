#Get taxa data overlap

dataset <- read.table("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/c70p3noAmbigs.snps", row.names=1, colClasses="character", skip=1)
dataset <- read.table("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/c5p3noAmbigs.snps", row.names=1, colClasses="character", skip=1)

source("~/phrynomics/trunk/phrynomicsFunctions.R")

#DataOverlap and MakePrresentAbsent function are stand alone in phrynomics. 

DataOverlap(dataset)[[2]]



#calculate branch length error from root to tips
#use BL.AllTrees
BL.AllTrees <- BL.AllTrees[[1]]

# GetAncestors, CalculateTotalTipBLError, and GetJustTipBLError function part of phrynomicsFunctions




CalculateTotalTipBLError(BL.AllTrees)

GetJustTipBLError(BL.AllTrees)






















