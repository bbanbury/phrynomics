################################################
##                                            ##
##          Phrynomics Functions              ##
##           edited: 6 May 2014               ##
##                                            ##
################################################

install.packages("~/phrynomics/phrynomics_1.0.tar.gz", type="source")
library(phrynomics)
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

files <- system("ls c*noAmbigs.snps", intern=T)
systemSeeds <- NULL
treeSeeds <- NULL
filenames <- NULL
pathToRAxML <- "raxmlHPC-PTHREADS-SSE3"
NumCores <- 16
models <- c("ASC_GTRCAT", "GTRCAT")
jobArgs <- NULL
for(model in sequence(length(models))) {
  for(i in sequence(length(files))) {
    inputFile <- files[i]
    outgroup <- "GAWI1"
    systemSeed <- floor(runif(1, min=1, max=10^6))
    bootstrapReps <- "autoMRE"
    treeSeed <- floor(runif(1, min=1, max=10^6))
    outputName <- paste(gsub("\\D+$", "", files[i]), "_", models[model], sep="")
    systemCall <- paste(pathToRAxML, " -T ", NumCores, " -s ", inputFile, " -m ", models[model], " -V -o ", outgroup, " -b ", systemSeed, " -# ", bootstrapReps, " -p ", treeSeed, " -n ", outputName, sep="")
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


#RAxMLfiles <- "~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/"
#MrBayesfiles <- "~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma/"

# Load RAxML Trees

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


# Load MrBayes Trees

analysis <- "MrBayes"
#setwd(MrBayesfiles)
trees <- system(paste("ls ", MrBayesfiles, "*.con.tre", sep=""), intern=T)  
TreeList <- CreateTreeList(trees, "MrBayes")
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


################################################
##     End Post-Analyses Tables and Figs      ##
################################################





































