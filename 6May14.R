##  ----------------------------------------  ##
##                                            ##
##          Phrynomics Functions              ##
##           edited: 6 May 2014               ##
##                                            ##
##  ----------------------------------------  ##

devtools::install_github("bbanbury/phrynomics")
library(phrynomics)
library(ape)

# Set location to where original pyrad files are located.
# This script will create 

setwd("/Users/Barb/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper")  



##  ----------------------------------------  ##
##            Begin File Creation             ##
##  ----------------------------------------  ##

#  Remove invariant sites and write phylip formatted data to file

snpFiles <- system("ls c*p3.snps", intern=T)
minSamples <- NULL
NuSNPs <- NULL
NuLoci <- NULL
for(i in sequence(length(snpFiles))) {
  filename <- strsplit(strsplit(snpFiles[i], "/")[[1]][length(strsplit(snpFiles[i], "/")[[1]])], ".", fixed=T)[[1]][1]
  snpList <- read.table(snpFiles[i], row.names=1, colClasses="character")
  snpList <- PreProcess(snpList)
  loci <- dim(snpList)[2]
  initialLociLengths <- nchar(snpList[1,])
  splitdata <- SplitSNP(snpList)
  KeepVector <- apply(splitdata, 2, IsVariable)  #If True, then keep as variable
  breaks <- which(splitdata[1,] == " ")
  write.table(rbind(splitdata, KeepVector), file=paste("KeepVector", filename, ".txt", sep=""))
  newSNP <- cSNP(splitdata, KeepVector=KeepVector)
  newLoci <- length(which(newSNP[1,] != ""))  #number of new loci 
  SNPFile <- paste(dim(newSNP)[1], length(which(KeepVector == "TRUE")) - length(breaks))
  write(SNPFile, file=paste(filename, "noAmbigs.snps", sep=""))
  write.table(newSNP, file=paste(filename, "noAmbigs.snps", sep=""), append=TRUE, quote=FALSE, col.names=FALSE)  
  minSamples <- c(minSamples, as.numeric(strsplit(filename, "\\D+", perl=T)[[1]][2]))
  NuSNPs <- c(NuSNPs, length(which(KeepVector == "TRUE")) - length(breaks))
  NuLoci <- c(NuLoci, newLoci)
}


# Writing RAxML invocation calls with bootstrapping
####change this to be exactly what the final script will write

files <- system("ls c*noAmbigs.snps", intern=T)
systemSeeds <- NULL
treeSeeds <- NULL
filenames <- NULL
pathToRAxML <- "raxmlHPC-PTHREADS-SSE3"
NumCores <- 6
models <- c("ASC_GTRCAT", "GTRCAT")
jobArgs <- NULL
for(model in sequence(length(models))) {
  for(i in sequence(length(files))) {
    inputFile <- files[i]
    outgroup <- "GAWI1"
    systemSeed <- floor(runif(1, min=1, max=10^6))
    bootstrapReps <- "50"
    treeSeed <- floor(runif(1, min=1, max=10^6))
    outputName <- paste(gsub("\\D+$", "", files[i]), "_", models[model], sep="")
    systemCall <- paste(pathToRAxML, " -T ", NumCores, " -s ", inputFile, " -f a -m ", models[model], " -V -o ", outgroup, " -# ", bootstrapReps, " -x ", systemSeed, " -p ", treeSeed, " -n ", outputName, sep="")
    jobArgs <- append(jobArgs, systemCall)
    systemSeeds <- c(systemSeeds, systemSeed) 
    treeSeeds <- c(treeSeeds, treeSeed) 
    filenames <- c(filenames, outputName) 
  }
save(systemSeeds, treeSeeds, filenames, file="seeds.rdata")
write(jobArgs, file=paste("RAXML_JobArgs", sep=""))
}


# Writing RAxML invocation calls for repeat runs (no bootstrapping)

files <- system("ls c*noAmbigs.snps", intern=T)
writejobArgs <- function(files, pathToRAxML="raxmlHPC-PTHREADS-SSE3", NumCores=6, NumRuns=20){
  models <- c("ASC_GTRCAT", "GTRCAT")
  jobArgs <- NULL
  for(whichFile in sequence(length(files))){ 
    file <- files[whichFile]
    for(model in sequence(length(models))){
      for(i in sequence(NumRuns)){
        outgroup <- "GAWI1"
        systemSeed <- floor(runif(1, min=1, max=10^6))
        treeSeed <- floor(runif(1, min=1, max=10^6))
        bootstrapReps <- 1
        outputName <- paste(gsub("\\D+$", "", file), "_", models[model], "_REP", i, sep="")
        systemCall <- paste(pathToRAxML, " -T ", NumCores, " -s ", file, " -m ", models[model], " -V -o ", outgroup, " -p ", treeSeed, " -n ", outputName, sep="")
        jobArgs <- append(jobArgs, systemCall)            
      }
    }
  }
  write(jobArgs, file=paste("repeat.RAXML_JobArgs", sep=""))
}
writejobArgs(files)


# Write MrBayes nexus formatted input files

files <- system("ls c*noAmbigs.snps", intern=T)
seedNumbers <- NULL
AllFilenames <- NULL
type <- c("variable", "all")
models <- c("ASC", "GTR")
for(model in sequence(length(models))){
  for(i in sequence(length(files))) {
    filename <- strsplit(strsplit(files[i], "/")[[1]][length(strsplit(files[i], "/")[[1]])], ".", fixed=T)[[1]][1]
    randomSeed <- floor(runif(1, min=1, max=100000))
    seedNumbers <- c(seedNumbers, randomSeed)
    snpList <- read.table(files[i], row.names=1, colClasses="character", skip=1)
    loci <- nchar(paste(snpList[1,], collapse=""))  #calculate number of loci BEFORE translating bases or the number of loci will be over estimated.  
    snpList2 <- TranslateBases(snpList, TRUE, "?")
    f <- paste(models[model], "_", filename, ".nex", sep="")  
    AllFilenames <- c(AllFilenames, f)
    write(paste("#NEXUS"), file=f)
    write(paste("[Written ", Sys.Date(), " via Rscript]", sep=""), file=f, append=TRUE)
    write(paste("BEGIN Data;"), file=f, append=TRUE)
    write(paste("   DIMENSIONS NTAX=", dim(snpList2)[1], " NCHAR=", loci, ";", sep=""), file=f, append=TRUE)
    write(paste("   FORMAT DATATYPE=Standard INTERLEAVE=no missing=?;", sep=""), file=f, append=TRUE)
    write(paste("Matrix"), file=f, append=TRUE)
    write.table(snpList2, file=f, append=TRUE, quote=FALSE, col.names=FALSE)  
    write(paste(""), file=f, append=TRUE)
    write(paste(";"), file=f, append=TRUE)
    write(paste("END;"), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("BEGIN mrbayes;"), file=f, append=TRUE)
    write(paste("outgroup GAWI1;"), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("set autoclose=yes;"), file=f, append=TRUE)
    write(paste("set seed=", randomSeed, ";", sep=""), file=f, append=TRUE)
    write(paste("lset rates=equal coding=", type[model], ";", sep=""), file=f, append=TRUE)
    write(paste("prset brlenspr=unconstrained:gammadir(1,1,1,1);"), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("mcmc nchains=4 ngen=20000000 printfreq=200000 samplefreq=10000;"), file=f, append=TRUE)
    write(paste("sumt conformat=simple burnin=500;", sep=""), file=f, append=TRUE)
    write(paste("sump burnin=500;", sep=""), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("END;"), file=f, append=TRUE)
  }
  save(seedNumbers, AllFilenames, file="seeds.Rdata")
}


##  ----------------------------------------  ##
##            End File Creation               ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##            Begin Run Analyses              ##
##  ----------------------------------------  ##

# Run RAxML via R (not parallelized)

toRun <- scan(file="RAXML_JobArgs", what="character", sep="\n")
for(i in sequence(length(toRun))){
  system(toRun[i])
}

toRun <- scan(file="repeat.RAXML_JobArgs", what="character", sep="\n")
for(i in sequence(length(toRun))){
  system(toRun[i])
}

# Run MrBayes via R (not parallelized)

toRun <- system("ls *.nex", intern=TRUE)
for(i in sequence(length(toRun))){
  system(paste("mb ", toRun[i], " > log", toRun[i], ".txt", sep=""))
}


##  ----------------------------------------  ##
##            End Run Analyses                ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##       Post-Analyses Trees and Data         ##
##  ----------------------------------------  ##


# Load RAxML Trees and post-analyses scraping

analysis <- "RAxML"
RAxML.trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstraps
RAxML.TreeList <- CreateTreeList(RAxML.trees, "RAxML")
treeMatrix <- CreateTreeMatrix(RAxML.trees)
treeMatrix2 <- AddTreeDist(treeMatrix, RAxML.TreeList)
RAxML.treeMatrix <- AddBLD(treeMatrix2, RAxML.TreeList)
BL.AllTrees.RAxML <- list()
for(i in sequence(dim(treeMatrix)[1])) {
  tree1 <- assTrees(RAxML.treeMatrix[i,1], RAxML.TreeList)[[1]]
  tree2 <- assTrees(RAxML.treeMatrix[i,2], RAxML.TreeList)[[1]]
  BL.AllTrees.RAxML[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  names(BL.AllTrees.RAxML)[[i]] <- rownames(treeMatrix)[i]
}
ML.results <- GetRAxMLStatsPostAnalysis(".")


# Load MrBayes Trees and post-analyses scraping

analysis <- "MrBayes"
MrBayes.trees <- system(paste("ls *.con.tre", sep=""), intern=T)  
MrBayes.TreeList <- CreateTreeList(MrBayes.trees, "MrBayes")
TreeList <- lapply(MrBayes.TreeList, multi2di)  #remove later
treeMatrix <- CreateTreeMatrix(MrBayes.trees)
treeMatrix2 <- AddTreeDist(treeMatrix, MrBayes.TreeList)
MrBayes.treeMatrix <- AddBLD(treeMatrix2, MrBayes.TreeList)
BL.AllTrees.MrBayes <- list()
for(i in sequence(dim(treeMatrix)[1])) {
  tree1 <- assTrees(MrBayes.treeMatrix[i,1], MrBayes.TreeList)[[1]]
  tree2 <- assTrees(MrBayes.treeMatrix[i,2], MrBayes.TreeList)[[1]]
  BL.AllTrees.MrBayes[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  names(BL.AllTrees.MrBayes)[[i]] <- rownames(treeMatrix)[i]
}
MB.results <- GetMrBayesStatsPostAnalysis(".")


# Create global objects for tables and figs to be made

levels <-NULL
for(i in sequence(length(files))){
  levels <- as.numeric(c(levels, strsplit(files[i], "\\D+")[[1]][2])) #numerical datasets from file names
}
orderedLevels <- sort(levels)  #numerical datasets from file names
whichDatasets <- paste("c", levels, "p3", sep="")  #datasets from file names
AllOrder <- paste("c", seq(5, 70, 5), "p3", sep="") #datasets from sequence
orderToGo <- AllOrder[AllOrder %in% whichDatasets]  #datasets of sequence that exist
focalDatasets <- c("c5p3", "c25p3", "c45p3", "c65p3")
whichFocalDatasets <- focalDatasets[focalDatasets %in% whichDatasets]
treeMatrices <- list(RAxML.treeMatrix, MrBayes.treeMatrix)
analyses <- c("RAxML", "MrBayes")


##  ----------------------------------------  ##
##      End Post-Analyses Trees and Data      ##
##  ----------------------------------------  ##






##  ----------------------------------------  ##
##            Make Tables and Figs            ##
##  ----------------------------------------  ##

# Table 1 was made by hand


# Make Table 2. Summary ddRadSeq data. 

dataset <- seq(from=70, to=5, by=-5)
ASC.ML.results <- ML.results[which(ML.results[,2] == "ASC"),]
table2 <- matrix(nrow=length(dataset), ncol=6)
rownames(table2) <- paste("c", dataset, "p3", sep="")
colnames(table2) <- c("Matrix", "Loci", "VariableSites", "Missing(%)", "PhrynoOverlap", "non-PhrynoOverlap")
table2[,1] <- paste("s", dataset, sep="")
rows <- match(ASC.ML.results[,1], rownames(table2))
table2 <- table2[!is.na(rows),]  
table2[,c(2,3,4)] <- as.matrix(ASC.ML.results[order(rows[which(!is.na(rows))]), c(3,4,6)])
dataOverlap <- list()
for(i in sequence(length(files))){
  dataset <- read.table(files[[i]], row.names=1, colClasses="character", skip=1)
  dataOverlap[[i]] <- DataOverlap(dataset)[[2]]
  names(dataOverlap)[[i]] <- strsplit(strsplit(files[i], "/")[[1]][length(strsplit(files[i], "/")[[1]])], "no")[[1]][1]
}
for(m in sequence(dim(table2)[1])){
  whichDataset <- which(names(dataOverlap) == rownames(table2)[m])
  table2[m,6] <- mean(dataOverlap[[whichDataset]][-grep(pattern="PH", names(dataOverlap[[whichDataset]]))])
  table2[m,5] <- mean(dataOverlap[[whichDataset]][grep(pattern="PH", names(dataOverlap[[whichDataset]]))])
}
write.table(table2, file="table2.txt", quote=FALSE, row.names=FALSE)  


# Table 3 was made by hand


# Figure 1 was made by hand


# Make Figure 2. Acquisition bias and tree length

pdf(file="Figure2.pdf", width=8.5, height=5)
layout(matrix(1:2, nrow=1, byrow=TRUE), respect=TRUE)
plot(c(orderedLevels, orderedLevels), ML.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data")
title(main="RAxML")
legend("topright", legend=c("GTR", "GTR + ASC"), col=c("blue", "lightblue"), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(orderToGo)-1)){
  dataToUse <- which(orderToGo[i] == ML.results$Level)
  nextDataToUse <- which(orderToGo[i+1] == ML.results$Level)
  segments(orderedLevels[i], ML.results$TreeLength[dataToUse[1]], orderedLevels[i+1], ML.results$TreeLength[nextDataToUse[1]], col="lightblue")
  segments(orderedLevels[i], ML.results$TreeLength[dataToUse[2]], orderedLevels[i+1], ML.results$TreeLength[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(orderToGo))){
  dataToUse <- which(orderToGo[i] == ML.results$Level)
  points(orderedLevels[i], ML.results$TreeLength[dataToUse[1]], pch=21, bg="lightblue")
  points(orderedLevels[i], ML.results$TreeLength[dataToUse[2]], pch=21, bg="blue")
}
Ymin <- min(MB.results$TreeLength.lowCI)
Ymax <- max(MB.results$TreeLength.uppCI)
plot(c(orderedLevels, orderedLevels), MB.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data", ylim=c(Ymin, Ymax))
title(main="MrBayes")
legend("topright", legend=c("GTR", "GTR + ASC"), col=c("blue", "lightblue"), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(orderToGo)-1)){
  dataToUse <- which(orderToGo[i] == MB.results$Level)
  nextDataToUse <- which(orderToGo[i+1] == MB.results$Level)
  segments(orderedLevels[i], MB.results$TreeLength[dataToUse[1]], orderedLevels[i+1], MB.results$TreeLength[nextDataToUse[1]], col="lightblue")
  segments(orderedLevels[i], MB.results$TreeLength[dataToUse[2]], orderedLevels[i+1], MB.results$TreeLength[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(orderToGo))){
  dataToUse <- which(orderToGo[i] == MB.results$Level)
  arrows(orderedLevels[i], MB.results$TreeLength.lowCI[dataToUse[1]], orderedLevels[i], MB.results$TreeLength.uppCI[dataToUse[1]], code=3, length=0.05, col="lightblue", angle=90)
  arrows(orderedLevels[i], MB.results$TreeLength.lowCI[dataToUse[2]], orderedLevels[i], MB.results$TreeLength.uppCI[dataToUse[2]], code=3, length=0.05, col="blue", angle=90)
  points(orderedLevels[i], MB.results$TreeLength[dataToUse[1]], pch=21, bg="lightblue")
  points(orderedLevels[i], MB.results$TreeLength[dataToUse[2]], pch=21, bg="blue")
}
dev.off()


# Make Figure 3. Scatterplot branch lengths

usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
pdf(file="Figure3.pdf", width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:8, nrow=2, byrow=TRUE), respect=TRUE)  
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  if(whichAnalysis == "RAxML")
    BL.AllTrees <- BL.AllTrees.RAxML
  if(whichAnalysis == "MrBayes")
    BL.AllTrees <- BL.AllTrees.MrBayes
  for(i in sequence(length(whichFocalDatasets))){
    dataToUse <- which(whichFocalDatasets[i] == names(BL.AllTrees))
    BLs <- BL.AllTrees[[dataToUse]]$branchlength[which(BL.AllTrees[[dataToUse]]$present)]
    corr.BLs <- BL.AllTrees[[dataToUse]]$corr.BL[which(BL.AllTrees[[dataToUse]]$present)]
    plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="non-ASC", type="n")
    linmod <- lm(corr.BLs ~ BLs)
    abline(linmod, lty=2)
    y <- 0.18
    points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
    points(BLs, corr.BLs, pch=21, bg="gray")
    text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
    lines(c(-1,1), c(-1,1))
    title(main=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""))
  }
}
dev.off()


# Figure 4 parts. Colored branch SNP phylogenies

figIndex <- 0
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  if(whichAnalysis == "RAxML"){
    TreeList <- RAxML.TreeList
    BL.AllTrees <- BL.AllTrees.RAxML
  }
  if(whichAnalysis == "MrBayes"){
    TreeList <- MrBayes.TreeList
    BL.AllTrees <- BL.AllTrees.MrBayes
  }
  for(i in sequence(length(orderToGo))){
    pdf(file=paste(whichAnalysis, ".", orderToGo[i], "trees.pdf", sep=""), width=8.5, height=11)
    figIndex <- figIndex+1
    dataToUse <- which(rownames(treeMatrices[[anal]]) == orderToGo[i])
    tree1 <- assTrees(treeMatrices[[anal]][dataToUse,1], TreeList)[[1]]
    tree1$tip.label[which(tree1$tip.label == "UMNO1")] <- "CADR2"  #change taxon
    tree1$tip.label[which(tree1$tip.label == "PHCO1")] <- "PHCE5"  #change taxon
    tree1$tip.label[which(tree1$tip.label == "PHCO3")] <- "PHCE6"  #change taxon
    tree2 <- assTrees(treeMatrices[[anal]][dataToUse,2], TreeList)[[1]]
    tree2$tip.label[which(tree2$tip.label == "UMNO1")] <- "CADR2" #change taxon in tree2
    tree2$tip.label[which(tree2$tip.label == "PHCO1")] <- "PHCE5" #change taxon in tree2
    tree2$tip.label[which(tree2$tip.label == "PHCO3")] <- "PHCE6" #change taxon in tree2
    plot(tree1, edge.lty=BL.AllTrees[[dataToUse]]$edgelty, edge.color=BL.AllTrees[[dataToUse]]$edgeColor, cex=0.5, edge.width=2)
    legtxt <- c("Discordant", "< -10%", "-10% to 10%", "> 10%", "> 20%", "> 30%", "> 40%", "> 50%")
    legcolors <- c("gray", rgb(51,51,255, max=255), "gray", rgb(255,255,102, max=255), rgb(255,178,102, max=255), rgb(225,128,0, max=255), rgb(225,0,0, max=255), rgb(153,0,0, max=255))
    legend("bottomleft", legend=legtxt, col=legcolors, lty=c(2,rep(1,7)), lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Branchlength Difference"))) 
    nodelabels(text=BL.AllTrees[[dataToUse]]$support[which(BL.AllTrees[[dataToUse]]$class == "internal")], node=BL.AllTrees[[dataToUse]]$desc[which(BL.AllTrees[[dataToUse]]$class == "internal")], cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, -0.1))
    nodelabels(text=BL.AllTrees[[dataToUse]]$corr.support[which(BL.AllTrees[[dataToUse]]$class == "internal")], node=BL.AllTrees[[dataToUse]]$desc[which(BL.AllTrees[[dataToUse]]$class == "internal")], cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, 1.1))
    if(whichAnalysis == "RAxML"){
      bquote1 <- bquote(bold("Supplemental Figure S" * .(as.character(figIndex)) * ".") * " Maximum likelihood phylogeny for the ascertainment-corrected dataset s" * .(strsplit(orderToGo[i], "\\D")[[1]][2]) ~ "(" * .(makeNumberwithcommas(ML.results[which(ML.results$Level == orderToGo[i]), 4][1])) ~ "variable sites," ~ "\n" )
      bquote2 <- bquote(.(unique(ML.results[which(ML.results$Level == orderToGo[i]), 6])) * "% missing data). Bootstrap values are shown on nodes (ascertainment-corrected above, no correction below). Branch")
      bquote3 <- bquote("color coding reflects the degree of relative length difference between ascertainment and non-ascertainment estimates.")
    }
    if(whichAnalysis == "MrBayes"){
      bquote1 <- bquote(bold("Supplemental Figure S" * .(as.character(figIndex)) * ".") * " Bayesian consensus phylogeny for the ascertainment-corrected dataset s" * .(strsplit(orderToGo[i], "\\D")[[1]][2]) ~ "(" * .(makeNumberwithcommas(MB.results[which(MB.results$Level == orderToGo[i]), 4][1])) ~ "variable sites," ~ "\n" )
      bquote2 <- bquote(.(unique(MB.results[which(MB.results$Level == orderToGo[i]), 6])) * "% missing data). Posterior probability values are shown on nodes (ascertainment-corrected above, no correction below). Branch")
      bquote3 <- bquote("color coding reflects the degree of relative length difference between ascertainment and non-ascertainment estimates.")
    }
    mtext(bquote1, side=3, cex=.75, adj=c(0), line=2)
    mtext(bquote2, side=3, cex=.75, adj=c(0), line=1)
    mtext(bquote3, side=3, cex=.75, adj=c(0), line=0)
    dev.off()
  }
}


# Make Figure 5. Mean branch length error (%)

pdf(file="Figure5.pdf", width=8.5, height=5)
for(anal in 1:length(analyses)){
  whichAnalysis <- analyses[anal] 
  if(whichAnalysis == "RAxML")
    BL.AllTrees <- BL.AllTrees.RAxML
  if(whichAnalysis == "MrBayes")
    BL.AllTrees <- BL.AllTrees.MrBayes
  taxon.BLdiff <- list()
  for(i in sequence(length(BL.AllTrees))){
    taxon.BLdiff[[i]] <- GetJustTipBLError(BL.AllTrees[[i]])
    names(taxon.BLdiff)[[i]] <- names(BL.AllTrees)[[i]]
  }
  mean.var.data <- NULL
  for(i in sequence(length(orderToGo))){
    whichData <- orderToGo[i]
    data1 <- taxon.BLdiff[[which(names(taxon.BLdiff) == whichData)]]
    data2 <- rep("OG", length(data1))
    data2[grep(pattern="PH", names(dataOverlap[[which(names(dataOverlap) == whichData)]]))] <- "PH"
    data3 <- data.frame(data1, data2)
    OGmean <- mean(data3[which(data3[,2] == "OG"),1])
    OGvar <- var(data3[which(data3[,2] == "OG"),1])
    PHmean <- mean(data3[which(data3[,2] == "PH"),1])
    PHvar <- var(data3[which(data3[,2] == "PH"),1])
    data4 <- c(whichData, OGmean, OGvar, PHmean, PHvar)
    mean.var.data <- rbind(mean.var.data, data4)
  }
  colnames(mean.var.data) <- c("whichData", "OGmean", "OGvar", "PHmean", "PHvar")
  mean.var.data <- as.data.frame(mean.var.data, row.names= mean.var.data[,1], stringsAsFactors=FALSE)
  for(i in 2:5){
    mean.var.data[,i] <- as.numeric(mean.var.data[,i])
  }
  if(any(rownames(mean.var.data) == "c70p3"))
    mean.var.data <- mean.var.data[-which(rownames(mean.var.data) == "c70p3"),]
  if(whichAnalysis == "RAxML")
    RAxML.mean.var.data <- mean.var.data
  if(whichAnalysis == "MrBayes")
    MrBayes.mean.var.data <- mean.var.data
}
q <- sequence(dim(mean.var.data)[1])
par(mar=c(5,5,2,5))
plot(rep(q, 2), c(RAxML.mean.var.data[,2], MrBayes.mean.var.data[,2]), type="n", ylim=c(min(c(MrBayes.mean.var.data[,4], MrBayes.mean.var.data[,2])), max(c(MrBayes.mean.var.data[,4], MrBayes.mean.var.data[,2]))), ylab="mean", xlab="", axes=FALSE)
axis(side=2)
axis(side=1, at=q, labels=orderToGo[-which(orderToGo == "c70p3")])
points(q, RAxML.mean.var.data[,2], col="gray")
points(q, RAxML.mean.var.data[,4], col="black")
points(q, MrBayes.mean.var.data[,2], col="gray")
points(q, MrBayes.mean.var.data[,4], col="black")
for(i in sequence(dim(RAxML.mean.var.data)[1]-1)){
  segments(i, RAxML.mean.var.data[i,2], i+1, RAxML.mean.var.data[i+1,2], col="gray", lty=1)
  segments(i, RAxML.mean.var.data[i,4], i+1, RAxML.mean.var.data[i+1,4], col="black", lty=1)
}
for(i in sequence(dim(MrBayes.mean.var.data)[1]-1)){
  segments(i, MrBayes.mean.var.data[i,2], i+1, MrBayes.mean.var.data[i+1,2], col="gray", lty=2)
  segments(i, MrBayes.mean.var.data[i,4], i+1, MrBayes.mean.var.data[i+1,4], col="black", lty=2)
}
dev.off()


# Make Figure 6 a-d. Scatterplot branch length support 
# Figure 6 e-h was done using AWTY online to account for values < 50% and bipartition files not matching up correctly. 

pdf(file="Figure6a-d.pdf", width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:4, nrow=1, byrow=TRUE), respect=TRUE)
whichAnalysis <- "RAxML"
treeMatrix <- RAxML.treeMatrix
BL.AllTrees <- BL.AllTrees.RAxML
for(i in sequence(length(focalDatasets))){
  dataToUse <- which(focalDatasets[i] == names(BL.AllTrees))
  datarows <- which(BL.AllTrees[[dataToUse]]$present)[which(BL.AllTrees[[dataToUse]]$present) %in% which(BL.AllTrees[[dataToUse]]$support != 0)]  #don't want all the tip support 0s or non homologous clades
  support <- BL.AllTrees[[dataToUse]]$support[datarows]
  corr.support <- BL.AllTrees[[dataToUse]]$corr.support[datarows]
  xlims <- ylims <- c(1,100)
  plot(support, corr.support, ylab="GTR Tree", xlab="ASC Tree", pch=21, bg="gray", xlim=xlims, ylim=ylims)
  
  #text(support, corr.support, labels=BL.AllTrees[[dataToUse]][datarows,1])
  title(main=paste("s", strsplit(focalDatasets[[i]], "\\D")[[1]][2], sep=""))
}
dev.off()


# Make Figure 7 a-b. Average RF distances among replicate runs.
# Figure 7c was made using Mesquite. 

MLRFdistMatrix <- GetRFmatrix("RAxML")
BIRFdistMatrix <- GetRFmatrix("MrBayes")
RFdistMatrix <- list(ML=MLRFdistMatrix, BI=BIRFdistMatrix)
pdf(file="Figure7a-b.pdf", width=5, height=8.5)
layout(matrix(1:2, nrow=2, byrow=TRUE), respect=TRUE)
titles <- c("Ave. RF - 20 Independent Runs", "Ave. RF - 20 Posterior Trees")
for(dist in sequence(length(RFdistMatrix))){
  plot(as.numeric(RFdistMatrix[[dist]][,3]), as.numeric(RFdistMatrix[[dist]][,4]), xlab="dataset", ylab="Average RF", xaxt="n", type="n", ylim=c(0,130))
  title(main= titles[[dist]])
  cols <- c("lightblue", "blue")
  pchs <- c(15, 19)
  axis(side=1, at=1:14, labels=RFdistMatrix[[dist]][1:14, 2])
  ascsub <- RFdistMatrix[[dist]][which(RFdistMatrix[[dist]][,1] == "ASC"),]
  gtrsub <- RFdistMatrix[[dist]][which(RFdistMatrix[[dist]][,1] == "GTR"),]
  points(as.numeric(ascsub[,3]), as.numeric(ascsub[,4]), col=cols[1], pch=pchs[[dist]])
  points(as.numeric(gtrsub[,3]), as.numeric(gtrsub[,4]), col=cols[2], pch=pchs[[dist]])
  for(i in 1:13){
    segments(as.numeric(ascsub[i,3]), as.numeric(ascsub[i,4]), as.numeric(ascsub[i+1,3]), as.numeric(ascsub[i+1,4]), col=cols[1])
    segments(as.numeric(gtrsub[i,3]), as.numeric(gtrsub[i,4]), as.numeric(gtrsub[i+1,3]), as.numeric(gtrsub[i+1,4]), col=cols[2])
  }
  legtxt <- c("ASC", "GTR")
  legcolors <- c(cols[1], cols[2])
  legend("topleft", legend=legtxt, col=legcolors, lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Model"))) 
}
dev.off()


# Figure 8 was made by hand



##  ----------------------------------------  ##
##         End Make Tables and Figs           ##
##  ----------------------------------------  ##





##  ----------------------------------------  ##
##      Double check RAxML Invocations        ##
##  ----------------------------------------  ##

if(analysis == "RAxML"){
  systemCallSeeds <- CheckInvocations(getwd())
  CheckSeeds(systemCallSeeds)
}



##  ----------------------------------------  ##
##     End Double check RAxML Invocations     ##
##  ----------------------------------------  ##
















