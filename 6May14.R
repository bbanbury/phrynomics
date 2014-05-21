################################################
##                                            ##
##          Phrynomics Functions              ##
##           edited: 6 May 2014               ##
##                                            ##
################################################

install.packages("~/phrynomics/phrynomics_1.2.tar.gz", type="source")
library(phrynomics)
library(ape)
setwd("~/Dropbox/UWstuff/phrynomics/")  #original pyRAD files located here



################################################
##            Begin File Creation             ##
################################################

#  Remove invariant sites and write phylip formatted data to file

snpFiles <- system("ls c*", intern=T)
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


# Writing RAxML invocation calls   
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


################################################
##            End File Creation               ##
################################################






################################################
##            Begin Run Analyses              ##
################################################

# Run RAxML via R (not parallelized)

toRun <- scan(file="RAXML_JobArgs", what="character", sep="\n")
for(i in sequence(length(toRun))){
  system(toRun[i])
}


# Run MrBayes via R (not parallelized)

toRun <- system("ls *.nex", intern=TRUE)
for(i in sequence(length(toRun))){
  system(paste("mb ", toRun[i], " > log", toRun[i], ".txt", sep=""))
}


################################################
##            End Run Analyses                ##
################################################






################################################
##      Post-Analyses Tables and Figs         ##
################################################


# Load RAxML Trees and post-analyses scraping

analysis <- "RAxML"
#setwd(RAxMLfiles)
trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstraps
TreeList <- CreateTreeList(trees, "RAxML")
treeMatrix <- CreateTreeMatrix(trees)
treeMatrix2 <- AddTreeDist(treeMatrix, TreeList)
treeMatrix3 <- AddBLD(treeMatrix2, TreeList)
BL.AllTrees.RAxML <- list()
for(i in sequence(dim(treeMatrix)[1])) {
  tree1 <- assTrees(treeMatrix3[i,1], TreeList)[[1]]
  tree2 <- assTrees(treeMatrix3[i,2], TreeList)[[1]]
  BL.AllTrees.RAxML[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  names(BL.AllTrees.RAxML)[[i]] <- rownames(treeMatrix)[i]
}
ML.results <- GetRAxMLStatsPostAnalysis(".")


# Load MrBayes Trees and post-analyses scraping

analysis <- "MrBayes"
#setwd(MrBayesfiles)
trees <- system(paste("ls *.con.tre", sep=""), intern=T)  
TreeList <- CreateTreeList(trees, "MrBayes")
TreeList <- lapply(TreeList, multi2di)  #remove later
treeMatrix <- CreateTreeMatrix(trees)
treeMatrix2 <- AddTreeDist(treeMatrix, TreeList)
treeMatrix3 <- AddBLD(treeMatrix2, TreeList)
BL.AllTrees.MrBayes <- list()
for(i in sequence(dim(treeMatrix)[1])) {
  tree1 <- assTrees(treeMatrix3[i,1], TreeList)[[1]]
  tree2 <- assTrees(treeMatrix3[i,2], TreeList)[[1]]
  BL.AllTrees.MrBayes[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  names(BL.AllTrees.MrBayes)[[i]] <- rownames(treeMatrix)[i]
}
MB.results <- GetMrBayesStatsPostAnalysis(".")


# Complete runs for indexing
levels <-NULL
for(i in sequence(length(files))){
  levels <- as.numeric(c(levels, strsplit(files[i], "\\D+")[[1]][2]))
}
whichDatasets <- paste("c", levels, "p3", sep="")  #datasets from file names
AllOrder <- paste("c", seq(5, 70, 5), "p3", sep="") #datasets from sequence
orderToGo <- AllOrder[AllOrder %in% whichDatasets]  #datasets of sequence that exist


# Make Table 2. Summary ddRadSeq data. 

dataset <- seq(from=70, to=5, by=-5)
ASC.ML.results <- ML.results[which(ML.results[,2] == "ASC"),]
table2 <- matrix(nrow=length(dataset), ncol=6)
rownames(table2) <- paste("c", dataset, "p3", sep="")
colnames(table2) <- c("Matrix", "Loci", "VariableSites", "Missing(%)", "PhrynoOverlap", "non-PhrynoOverlap")
table2[,1] <- paste("s", dataset, sep="")
rows <- match(rownames(table2), ASC.ML.results[,1])
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


# Make Figure 2. Acquisition bias and tree length

layout(matrix(1:2, nrow=1, byrow=TRUE), respect=TRUE)
plot(c(levels, levels), ML.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data")
title(main="RAxML")
legend("topright", legend=c("GTR", "GTR + ASC"), col=c("blue", "lightblue"), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(whichDatasets)-1)){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  nextDataToUse <- which(whichDatasets[i+1] == ML.results$Level)
  segments(levels[i], ML.results$TreeLength[dataToUse[1]], levels[i+1], ML.results$TreeLength[nextDataToUse[1]], col="lightblue")
  segments(levels[i], ML.results$TreeLength[dataToUse[2]], levels[i+1], ML.results$TreeLength[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  points(levels[i], ML.results$TreeLength[dataToUse[1]], pch=21, bg="lightblue")
  points(levels[i], ML.results$TreeLength[dataToUse[2]], pch=21, bg="blue")
}
Ymin <- min(MB.results$TreeLength.lowCI)
Ymax <- max(MB.results$TreeLength.uppCI)
plot(c(levels, levels), MB.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data", ylim=c(Ymin, Ymax))
title(main="MrBayes")
legend("topright", legend=c("GTR", "GTR + ASC"), col=c("blue", "lightblue"), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(whichDatasets)-1)){
  dataToUse <- which(whichDatasets[i] == MB.results$Level)
  nextDataToUse <- which(whichDatasets[i+1] == MB.results$Level)
  segments(levels[i], MB.results$TreeLength[dataToUse[1]], levels[i+1], MB.results$TreeLength[nextDataToUse[1]], col="lightblue")
  segments(levels[i], MB.results$TreeLength[dataToUse[2]], levels[i+1], MB.results$TreeLength[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == MB.results$Level)
  arrows(levels[i], MB.results$TreeLength.lowCI[dataToUse[1]], levels[i], MB.results$TreeLength.uppCI[dataToUse[1]], code=3, length=0.05, col="lightblue", angle=90)
  arrows(levels[i], MB.results$TreeLength.lowCI[dataToUse[2]], levels[i], MB.results$TreeLength.uppCI[dataToUse[2]], code=3, length=0.05, col="blue", angle=90)
  points(levels[i], MB.results$TreeLength[dataToUse[1]], pch=21, bg="lightblue")
  points(levels[i], MB.results$TreeLength[dataToUse[2]], pch=21, bg="blue")
}


# Make Figure 3. Scatterplot branch lengths



















# Make Figure 4. Colored branch SNP phylogenies
##DOESN'T WORK YET!!!

for(analysis in c("RAxML", "MrBayes")){ 
  for(i in sequence(length(orderToGo))){
    letters <- as.character(sequence(length(orderToGo)))
    dataToUse <- which(rownames(treeMatrix) == orderToGo[i])
    pdf(file=paste(analysis, ".", rownames(treeMatrix)[dataToUse], "trees.pdf", sep=""), width=8.5, height=11)
    tree1 <- assTrees(treeMatrix3[dataToUse,1], TreeList)[[1]]
    tree1$tip.label[which(tree1$tip.label == "UMNO1")] <- "CADR2"  #change taxon
    tree2 <- assTrees(treeMatrix3[dataToUse,2], TreeList)[[1]]
    tree2$tip.label[which(tree2$tip.label == "UMNO1")] <- "CADR2" #change taxon in tree2
    plot(tree1, edge.lty=BL.AllTrees[[dataToUse]]$edgelty, edge.color=BL.AllTrees[[dataToUse]]$edgeColor, cex=0.5, edge.width=2)
    legtxt <- c("Discordant", "< -10%", "-10% to 10%", "> 10%", "> 20%", "> 30%", "> 40%", "> 50%")
    legcolors <- c("gray", rgb(51,51,255, max=255), "gray", rgb(255,255,102, max=255), rgb(255,178,102, max=255), rgb(225,128,0, max=255), rgb(225,0,0, max=255), rgb(153,0,0, max=255))
    legend("bottomleft", legend=legtxt, col=legcolors, lty=c(2,rep(1,7)), lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Branchlength Difference"))) 
    nodelabels(text=BL.AllTrees[[dataToUse]]$support[which(BL.AllTrees[[dataToUse]]$class == "internal")], node=BL.AllTrees[[dataToUse]]$desc[which(BL.AllTrees[[dataToUse]]$class == "internal")], cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, -0.1))
    nodelabels(text=BL.AllTrees[[dataToUse]]$corr.support[which(BL.AllTrees[[dataToUse]]$class == "internal")], node=BL.AllTrees[[dataToUse]]$desc[which(BL.AllTrees[[dataToUse]]$class == "internal")], cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, 1.1))
    if(analysis == "RAxML"){
      bquote1 <- bquote(bold("Supplemental Figure S" * .(letters[i]) * ".") * " Maximum likelihood phylogeny for the ascertainment-corrected dataset s" * .(strsplit(orderToGo[i], "\\D")[[1]][2]) ~ "(" * .(makeNumberwithcommas(RAxMLresults[which(RAxMLresults$Level == orderToGo[i]), 4][1])) ~ "variable sites," ~ "\n" )
      bquote2 <- bquote(.(unique(RAxMLresults[which(RAxMLresults$Level == orderToGo[i]), 6])) * "% missing data). Bootstrap values are shown on nodes (ascertainment-corrected above, no correction below). Branch")
      bquote3 <- bquote("color coding reflects the degree of relative length difference between ascertainment and non-ascertainment estimates.")
    }
    if(analysis == "MrBayes"){
      bquote1 <- bquote(bold("Supplemental Figure S2" * .(letters[i]) * ".") * " Bayesian consensus phylogeny for the ascertainment-corrected dataset s" * .(strsplit(orderToGo[i], "\\D")[[1]][2]) ~ "(" * .(makeNumberwithcommas(RAxMLresults[which(RAxMLresults$Level == orderToGo[i]), 4][1])) ~ "variable sites," ~ "\n" )
      bquote2 <- bquote(.(unique(RAxMLresults[which(RAxMLresults$Level == orderToGo[i]), 6])) * "% missing data). Posterior probability values are shown on nodes (ascertainment-corrected above, no correction below). Branch")
      bquote3 <- bquote("color coding reflects the degree of relative length difference between ascertainment and non-ascertainment estimates.")
    }
    mtext(bquote1, side=3, cex=.75, adj=c(0), line=2)
    mtext(bquote2, side=3, cex=.75, adj=c(0), line=1)
    mtext(bquote3, side=3, cex=.75, adj=c(0), line=0)
    dev.off()
  }
}








# Make Figure 5. Mean branch length error (%)

# Make Figure 6. Scatterplot branch length support 





################################################
##     End Post-Analyses Tables and Figs      ##
################################################

