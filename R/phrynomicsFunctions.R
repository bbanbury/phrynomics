Get##Set of functions to deal with SNP data and compare results
##Barb Banbury  14Nov13
##RAxMLsetup and RAxMLanalysis for the SNP setup and tree analysis

library(ape)
library(phangorn)


################################################
## Code to deal with SNP data in various ways ##
################################################

#code removed 4-23 and split off to R/ 


################################################
#################   END   ######################
################################################



################################################
####  Code to deal with data/tree analysis #####
################################################

getMissingDataAmount <- function(filename){
  split1 <- strsplit(filename, "\\D+")[[1]][2]
  return(as.numeric(split1))
}

CreateTreeList <- function(filenames, analysis="RAxML"){
  analysis <- match.arg(arg=analysis, choices=c("RAxML", "MrBayes"), several.ok=FALSE)
  TreeList <- list()
  if(analysis == "RAxML"){
    for(i in sequence(length(filenames))){
      TreeList[[i]] <- ladderize(read.tree(filenames[i]))
      names(TreeList)[[i]] <- filenames[i]
    }
  }
  if(analysis == "MrBayes"){
    for(i in sequence(length(filenames))){
      TreeList[[i]] <- ladderize(read.nexus(filenames[i])[[1]])
      names(TreeList)[[i]] <- filenames[i]
    }
  }
  return(TreeList)
}

CreateTreeMatrix <- function(trees) {
#Creates a matrix of file names that correspons to missing data amounts (rows) and model (cols)
  missingDataTypes <- sapply(trees, getMissingDataAmount)
  ASCtrees <- grep("ASC", trees)
  GTRtrees <- grep("GTR", trees)
  treeMatrix <- matrix(nrow=length(unique(missingDataTypes)), ncol=2)
  rownames(treeMatrix) <- paste("c", unique(missingDataTypes), "p3", sep="")
  colnames(treeMatrix) <- c("ASC", "GTR")
  for(i in sequence(dim(treeMatrix)[1])) {
    twoFiles <- trees[grep(rownames(treeMatrix)[i], trees, fixed=TRUE)]
    treeMatrix[i,] <- c(twoFiles[grep("ASC", twoFiles)], twoFiles[-grep("ASC", twoFiles)])
  }
  return(as.data.frame(treeMatrix, stringsAsFactors=FALSE))
}

assTrees <- function(TreeMatrixName, ListOfTrees){
  return(ListOfTrees[which(names(ListOfTrees) == TreeMatrixName)])
}

AddTreeDist <- function(TreeMatrixName, ListOfTrees){
#Uses phangorn to calculate Steel & Penny whole-tree metrics
#Steel M. A. and Penny P. (1993) Distributions of tree comparison metrics - some new results, Syst. Biol.,42(2), 126-141
  treeDists <- NULL
  for(row in sequence(dim(TreeMatrixName)[1])){
    twoTrees <- c(assTrees(TreeMatrixName[row,1], ListOfTrees), assTrees(TreeMatrixName[row,2], ListOfTrees))
    dists <- phangorn::treedist(twoTrees[[1]], twoTrees[[2]])
    treeDists <- rbind(treeDists, dists)
  }
  cbind(TreeMatrixName, treeDists)
}

##  Kuhner & Felsenstein Branch Length distance measure (Kscore)
AddBLD <- function(TreeMatrixName, ListOfTrees){
#Kuhner, M. K. and J. Felsenstein. 1994. A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates. Molecular Biology and Evolution 11: 459-468.
#This uses Soria-Carrasco & Jose Castresana perl script to calculate K-score via terminal
#Soria-Carrasco, V., Talavera, G., Igea, J., and Castresana, J. (2007). The K tree score: quantification of differences in the relative branch length and topology of phylogenetic trees. Bioinformatics 23, 2954-2956.
#http://molevol.cmima.csic.es/castresana/Ktreedist.html
  Kscores <- NULL
  for (row in sequence(dim(TreeMatrixName)[1])) {
    twoTrees <- c(assTrees(TreeMatrixName[row,1], ListOfTrees), assTrees(TreeMatrixName[row,2], ListOfTrees))
    twoTrees[[1]]$node.label <- NULL
    twoTrees[[2]]$node.label <- NULL
    write.nexus(twoTrees[[1]], file="tree1")  #prog will not run with bootstrap vals
    write.nexus(twoTrees[[2]], file="tree2")  
    runProg <- system("perl /Applications/PhylogeneticsPrograms/Ktreedist_v1/Ktreedist.pl -rt tree1 -ct tree2 -a", intern=T)
    vals <- grep("UNTITLED", runProg)[2]
    kscore <- as.numeric(strsplit(runProg[vals], " ")[[1]][grep("\\d+", strsplit(runProg[vals], " ")[[1]])])
    names(kscore) <- c("Kscore", "ScaleFactor", "SymmDiff", "Npartitions")
    Kscores <- rbind(Kscores, kscore)
  }
  cbind(TreeMatrixName, Kscores)
}


GetEdgeList <- function(tree) {
  desc.order <- tree$edge[,2]  #store this order so that you can go back to it after
  tipList <- cbind(tree$edge, tree$edge[, 2] %in% tree$edge[, 1], tree$edge.length)
  tipList <- as.data.frame(tipList, stringsAsFactors=F)
  colnames(tipList) <- c("anc", "desc", "class", "branchlength")
  tipList[which(tipList[,3] == 0), 3] <- "tip"
  tipList[which(tipList[, 3] == 1), 3] <- "internal"
  tipList$support <- rep(0, dim(tipList)[1])
  tipList  <- tipList[order(tipList[,2]), ]  # reorder edge list so that you can assign the correct support vals
  if(!is.null(tree$node.label))
    tipList$support[which(tipList$class == "internal")] <- as.numeric(tree$node.label)[-1] #node support comes in with an extra space in the beginning, so it has to be cleaved then readded for plotting.
  tipList <- tipList[match(desc.order , tipList[,2]),]
  options(digits=15)
  tipList$branchlength <- as.numeric(tipList$branchlength)
  return(data.frame(tipList, stringsAsFactors = F))
}


nodeOffspring <- function(tree, anc.node) {
  rows <- which(tree$edge[, 1] == anc.node)
  return(tree$edge[rows, 2])
}

nodeLeaves <- function(tree, anc.node) {
  ntips <- Ntip(tree)
  if (anc.node <= ntips) 
    return(tree$tip.label[as.numeric(anc.node)])
  listTaxa <- character()
  descendents <- nodeOffspring(tree, anc.node)
  for (j in descendents) {
    if (j <= ntips) 
      listTaxa <- c(listTaxa, tree$tip.label[as.numeric(j)])
    else listTaxa <- c(listTaxa, nodeLeaves(tree, j))  #recursive....
  }
  return(listTaxa)
}

CheckSharedMonophy <- function(t1.node, tree1, tree2) {
#t1.node should be desc.nodes--a single edge (ex: 74)
#returns T/F if the taxa from that edge match a monophyletic group in tree 2  
  t1.taxa <- nodeLeaves(tree1, t1.node)
  if(length(t1.taxa) == 1)
    return(TRUE)
  mrca <- getMRCA(tree2, t1.taxa)
  t2.taxa <- nodeLeaves(tree2, mrca)
  if(all(t1.taxa %in% t2.taxa) && all(t2.taxa %in% t1.taxa))
    return(TRUE)
  else return(FALSE)  
}

GetCorrespondingEdge <- function(t1.node, tree1, tree2) {
  if(CheckSharedMonophy(t1.node, tree1, tree2)) {
    t1.taxa <- nodeLeaves(tree1, t1.node)
    if(length(t1.taxa) == 1)
      return(which(tree2$tip.label == t1.taxa))
    return(getMRCA(tree2, t1.taxa))
  }
  else return(0)
}
#GetCorrespondingEdge(77, tree1, tree2)

GetCorresonding <- function(corr.desc, t2){  
  return(t2[which(t2$desc == corr.desc),])
}

GetSingleEdgeColor <- function(relativeBLdiff) {
  if(is.na(relativeBLdiff))  return(NA)
  else if (relativeBLdiff < -10) return(rgb(51,51,255, maxColorValue =255)) #underestimate over 10%
  else if (relativeBLdiff <= 10) return("gray")  #plus/minus 10%
  else if (relativeBLdiff < 20) return(rgb(255,255,102, maxColorValue=255))
  else if (relativeBLdiff < 30) return(rgb(255,178,102, maxColorValue=255))
  else if (relativeBLdiff < 40) return(rgb(225,128,0, maxColorValue=255))
  else if (relativeBLdiff < 50) return(rgb(225,0,0, maxColorValue=255))
  else return(rgb(153,0,0, maxColorValue=255))	  
}
#GetSingleEdgeColor(30)

ReturnMinCI <- function(x, vstat){
  if(x != 0)
    return(vstat[which(vstat$Median == x),3])
  else return(0)
}

ReturnMaxCI <- function(x, vstat){
  if(x != 0)
    return(vstat[which(vstat$Median == x),4])
  else return(0)
}


MakeBranchLengthMatrix <- function(tree1, tree2, analysis="RAxML", dataset=NULL){  
  t1 <- GetEdgeList(tree1) 
  t2 <- GetEdgeList(tree2) 
  desc.nodes <- t1[,2]
  t1$present <- sapply(desc.nodes, CheckSharedMonophy, tree1=tree1, tree2=tree2) 
  t1$corr.anc <- sapply(t1[,1], GetCorrespondingEdge, tree1=tree1, tree2=tree2)
  t1$corr.desc <- sapply(desc.nodes, GetCorrespondingEdge, tree1=tree1, tree2=tree2)
  t1$corr.BL <- rep(0, dim(t1)[1])
  t1$corr.support <- rep(0, dim(t1)[1])
  for(i in which(t1$present)){
    t1$corr.BL[i] <- GetCorresonding(t1$corr.desc[i], t2)[[4]]
    t1$corr.support[i] <- GetCorresonding(t1$corr.desc[i], t2)[[5]]    
  }
  t1$BL.DIFF <- t1$corr.BL - t1$branchlength #tree B - tree A
  t1$BL.DIFF[which(!t1$present)] <- 0  #remove data when not comparable
  t1$relativeBLdiff <- (t1$BL.DIFF / t1$branchlength) * 100
  t1$edgeColor <- sapply(t1$relativeBLdiff, GetSingleEdgeColor)
  t1$edgelty <- rep(3, dim(t1)[1])
  t1$edgelty[which(t1$present)] <- 1
  if(analysis == "MrBayes"){
    vstat1 <- read.table(system(paste("ls ASC_", dataset, "*.vstat", sep=""), intern=TRUE), row.names=1, stringsAsFactors=FALSE, skip=1)
    vstat2 <- read.table(system(paste("ls GTR_", dataset, "*.vstat", sep=""), intern=TRUE), row.names=1, stringsAsFactors=FALSE, skip=1)
    colnames(vstat1) <- colnames(vstat2) <- read.table(system(paste("ls ASC_", dataset, "*.vstat", sep=""), intern=TRUE), row.names=1, stringsAsFactors=FALSE)[1,]
    t1$min.ASC.BL.CI <- sapply(t1$branchlength, ReturnMinCI, vstat=vstat1)
    t1$max.ASC.BL.CI <- sapply(t1$branchlength, ReturnMaxCI, vstat=vstat1)
    t1$min.GTR.BL.CI <- sapply(t1$corr.BL, ReturnMinCI, vstat=vstat2)
    t1$max.GTR.BL.CI <- sapply(t1$corr.BL, ReturnMaxCI, vstat=vstat2)
  }
  return(t1)
}

makeNumberwithcommas <- function(number){
  string <- strsplit(as.character(number), "")[[1]]
  if(length(string) > 3){
    string <- paste(paste(paste(string[1:(length(string)-3)], collapse=""), ",", sep=""), paste(string[-1:-(length(string)-3)], collapse=""), sep="")
  }
  else string <- paste(string, collapse="")
  return(string)
}

GetSupportValue <- function(lineNumber, tstat.table){
  lineNumber <- as.numeric(lineNumber)
  if(lineNumber %in% tstat.table[,1])
    return(tstat.table[which(tstat.table[,1] == lineNumber), 3])
}


CompareMrBayesPosteriors <- function(run1, run2, ntax=73){
#run1 and run2 should correspond with ASC and GTR models
#should be able to give it the treeMatrix and have it find the associated files
  ASCbipartsFile <- read.table(system(paste("ls ",strsplit(run1, ".", fixed=TRUE)[[1]][1], "*.parts", sep=""), intern=TRUE), stringsAsFactors=FALSE, skip=1)
  ASCbipartsFile <- ASCbipartsFile[-which(ASCbipartsFile == "ID"):-dim(ASCbipartsFile)[1],]
  GTRbipartsFile <- read.table(system(paste("ls ",strsplit(run2, ".", fixed=TRUE)[[1]][1], "*.parts", sep=""), intern=TRUE), skip=1, stringsAsFactors=FALSE)
  GTRbipartsFile <- GTRbipartsFile[-which(GTRbipartsFile == "ID"):-dim(GTRbipartsFile)[1],]
  ASCtstatFile <- read.table(system(paste("ls ",strsplit(run1, ".", fixed=TRUE)[[1]][1], "*.tstat", sep=""), intern=TRUE), skip=1, stringsAsFactors=FALSE)
  GTRtstatFile <- read.table(system(paste("ls ",strsplit(run2, ".", fixed=TRUE)[[1]][1], "*.tstat", sep=""), intern=TRUE), skip=1, stringsAsFactors=FALSE)  
  ASC.bipartInternalNodes <- (ntax+1):max(as.numeric(ASCbipartsFile[,1]))
  ASC.results <- matrix(nrow=length(ASC.bipartInternalNodes), ncol=5)
  ASC.results[,1] <- ASC.bipartInternalNodes
  ASC.results[,2] <- sapply(ASC.bipartInternalNodes, GetSupportValue, tstat.table=ASCtstatFile)
  for(line in sequence(length(ASC.bipartInternalNodes))){
    pattern <- ASCbipartsFile[which(ASCbipartsFile[,1] == ASC.bipartInternalNodes[line]),2]
    ASC.results[line,3] <- pattern
    corr.line <- GTRbipartsFile[which(GTRbipartsFile[,2] == pattern), 1]
    if(length(corr.line) == 0){
      ASC.results[line,4] <- NA
      ASC.results[line,5] <- 0
    }
    if(length(corr.line) > 0){
      ASC.results[line,4] <- corr.line
      ASC.results[line,5] <- GetSupportValue(corr.line, GTRtstatFile)
    }
  }
  GTR.bipartInternalNodes <- (ntax+1):max(as.numeric(GTRbipartsFile[,1]))
  GTR.results <- matrix(nrow=length(GTR.bipartInternalNodes), ncol=5)
  GTR.results[,4] <- GTR.bipartInternalNodes
  GTR.results[,5] <- sapply(GTR.bipartInternalNodes, GetSupportValue, tstat.table=GTRtstatFile)
  for(line in sequence(length(GTR.bipartInternalNodes))){
    pattern <- GTRbipartsFile[which(GTRbipartsFile[,1] == GTR.bipartInternalNodes[line]),2]
    GTR.results[line,3] <- pattern
    corr.line <- ASCbipartsFile[which(ASCbipartsFile[,2] == pattern), 1]
    if(length(corr.line) == 0){
      GTR.results[line,1] <- NA
      GTR.results[line,2] <- 0
    }
    if(length(corr.line) > 0){
      GTR.results[line,1] <- corr.line
      GTR.results[line,2] <- GetSupportValue(corr.line, ASCtstatFile)
    }
  }
results <- as.data.frame(rbind(ASC.results, GTR.results), stringsAsFactors=FALSE)
colnames(results) <- c("ASC.bipart.line", "ASC.support", "pattern", "GTR.bipart.line", "GTR.support")
results$ASC.support <- as.numeric(results$ASC.support)
results$GTR.support <- as.numeric(results$GTR.support)
return(results)
}

getColor <- function(BLtable, nonSigColor="gray", sigColor="red", method=c("varOverlap", "meanWithin")){
#for method, choose whether you want to distinguish significant deviants either by non-overlapping CI distributions or by the mean of the GTR BL not falling within the ASC distribution 
  method <- match.arg(method, choices=c("varOverlap", "meanWithin"))  
  if(is.null(BLtable$min.ASC.BL.CI))
    return(rep(nonSigColor, dim(BLtable)[1]))
  colorVector <- rep("NA", dim(BLtable)[1])
  for(i in sequence(dim(BLtable)[1])){
    if(BLtable$corr.BL[i] != 0){
      if(method == "meanWithin"){
        if(BLtable$min.ASC.BL.CI[i] < BLtable$corr.BL[i] && BLtable$corr.BL[i] < BLtable$max.ASC.BL.CI[i])
          colorVector[i] <- nonSigColor
        else colorVector[i] <- sigColor
      }
      if(method == "varOverlap"){
        if(BLtable$min.GTR.BL.CI[i] < BLtable$max.ASC.BL.CI[i] || BLtable$min.ASC.BL.CI[i] < BLtable$max.GTR.BL.CI[i])
          colorVector[i] <- nonSigColor
        else  colorVector[i] <- sigColor
      }
    }
  }
return(colorVector)
}

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



################################################
#################   END   ######################
################################################




################################################
#######    Post Analysis Scraping   ############
################################################

GetRAxMLStatsPostAnalysis <- function(workingDirectoryOfResults) {
  startingDir <- getwd()
  setwd(workingDirectoryOfResults)
  cFiles <- system("ls c*noAmbigs.snps", intern=T)
  outFiles <- system("ls RAxML_info*", intern=T)
  results <- matrix(nrow=length(outFiles), ncol=10)
  for(i in sequence(length(outFiles))){
    MissingDataLevel <- strsplit(strsplit(outFiles[i], "_")[[1]][2], "info.")[[1]][2]
    whichModel <- strsplit(outFiles[i], "_")[[1]][3]
    SNPdataset <- read.table(cFiles[grep(MissingDataLevel, cFiles)], row.names=1, colClasses="character", skip=1)
    VariableSites <- nchar(paste(SNPdataset[1,], collapse=""))
    numberLoci <- dim(SNPdataset)[2]
    alignmentPatterns <- gsub("\\D", "", system(paste("grep 'distinct alignment patterns'", outFiles[i]), intern=T))
    Missing <-gsub("[A-Za-z:]+|[%]$", "", system(paste("grep 'Proportion of gaps and completely undetermined characters in this alignment:'", outFiles[i]), intern=T), perl=T)
    BootstrapTime <- strsplit(system(paste("grep 'Overall Time for '", outFiles[i]), intern=T), split="[A-Za-z:]+|[%]$", perl=T)[[1]][6]
    Likelihood <- strsplit(system(paste("grep 'Final ML Optimization Likelihood:'", outFiles[i]), intern=T), split="[A-Za-z:]+|[%]$", perl=T)[[1]][5]
    #Alpha <- gsub("[alpha: ]", "", system(paste("grep alpha: ", outFiles[i]), intern=T), perl=T)
    Alpha <- "0"
	TreeLength <- gsub("Tree-Length: ", "", system(paste("grep Tree-Length: ", outFiles[i]), intern=T), fixed=T)
    results[i,] <- c(MissingDataLevel, whichModel, numberLoci, VariableSites, alignmentPatterns, Missing, BootstrapTime, Likelihood, Alpha, TreeLength)
  }
  results <- data.frame(results, stringsAsFactors=FALSE)
  colnames(results) <- c("Level", "Model", "NumberLoci", "VariableSites", "AlignmentPatterns", "MissingData", "BootstrapTime", "Likelihood", "Alpha", "TreeLength")
  options(digits=10)
  results$Level <- as.factor(results$Level)
  results$Model <- as.factor(results$Model)
  results$NumberLoci <- as.numeric(results$NumberLoci)
  results$VariableSites <- as.numeric(results$VariableSites)
  results$AlignmentPatterns <- as.numeric(results$AlignmentPatterns)
  results$MissingData <- as.numeric(results$MissingData)
  results$BootstrapTime <- as.numeric(results$BootstrapTime)
  results$Likelihood <- as.numeric(results$Likelihood)
  results$Alpha <- as.numeric(results$Alpha)
  results$TreeLength <- as.numeric(results$TreeLength)
  setwd(startingDir)
  return(results)
}

GetLinesToSkip <- function(file){
  return(length(system(paste("grep 'ID:' ", file, sep=""), intern=TRUE)))
}

GetMrBayesStatsPostAnalysis <- function(workingDirectoryOfResults){
  startingDir <- getwd()
  setwd(workingDirectoryOfResults)
  nexFiles <- system("ls *.nex", intern=TRUE)  #num loci
  pstatFiles <- system("ls *.pstat", intern=TRUE)  #tree length and alpha
  logFiles <- system("ls log*", intern=TRUE)  
  results <- matrix(nrow=length(nexFiles), ncol=12)
  for(i in sequence(length(nexFiles))){
    MissingDataLevel <- paste("c", strsplit(nexFiles[i], "\\D+")[[1]][2], "p3", sep="")
    whichModel <- strsplit(nexFiles[i], "_")[[1]][1]
    numLociLine <- system(paste("grep DIMENSIONS", nexFiles[i]), intern=TRUE)
    numberLoci <- strsplit(numLociLine, "\\D+")[[1]][3]
    dataName <- paste(whichModel, "_c", strsplit(nexFiles[i], "\\D+")[[1]][2], "p3", sep="")
    pstats <- read.csv(pstatFiles[grep(dataName, pstatFiles)], sep="", skip=GetLinesToSkip(pstatFiles[grep(dataName, pstatFiles)]))
    treelength <- pstats[1,2]
    treelength.lowCI <- pstats[1,4]
    treelength.uppCI <- pstats[1,5]
    treelengthESS <- pstats[1,8]
    alpha <- pstats[2,2]
    alpha.lowCI <- pstats[2,4]
    alpha.uppCI <- pstats[2,5]
    alphaESS <- pstats[2,8]
    stdSplitLine <- system(paste("grep -A 1 'Summary statistics for partitions with frequency' ", logFiles[grep(dataName, pstatFiles)], sep=""), intern=T)[2]

    stdSplits <- gsub("          Average standard deviation of split frequencies = ", "", stdSplitLine)
    results[i,] <- c(MissingDataLevel, whichModel, numberLoci, treelength, treelength.lowCI, treelength.uppCI, treelengthESS, alpha, alpha.lowCI, alpha.uppCI, alphaESS, stdSplits)
  }  
  results <- data.frame(results, stringsAsFactors=FALSE)
  colnames(results) <- c("Level", "Model", "NumberLoci", "TreeLength", "TreeLength.lowCI", "TreeLength.uppCI", "TreeLengthESS", "Alpha", "Alpha.lowCI", "Alpha.uppCI", "AlphaESS", "stdSplits")
  results$NumberLoci <- as.numeric(results$NumberLoci)
  results$TreeLength <- as.numeric(results$TreeLength)
  results$TreeLength.lowCI <- as.numeric(results$TreeLength.lowCI)
  results$TreeLength.uppCI <- as.numeric(results$TreeLength.uppCI)
  results$TreeLengthESS <- as.numeric(results$TreeLengthESS)
  results$Alpha <- as.numeric(results$Alpha)
  results$Alpha.lowCI <- as.numeric(results$Alpha.lowCI)
  results$Alpha.uppCI <- as.numeric(results$Alpha.uppCI)
  results$AlphaESS <- as.numeric(results$AlphaESS)
  results$stdSplits <- as.numeric(results$stdSplits)
  setwd(startingDir)
  return(results)
}





################################################
#################   END   ######################
################################################








