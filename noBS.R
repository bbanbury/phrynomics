
# Writing RAxML invocation calls for multiple ML searches, NO BS
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLrepeatMLsearch")
files <- system("ls c*.snps", intern=TRUE)
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
        outputName <- paste(gsub("\\D+$", "", file), "_", models[model], "_", i, sep="")
        systemCall <- paste(pathToRAxML, " -T ", NumCores, " -s ", file, " -m ", models[model], " -V -o ", outgroup, " -p ", treeSeed, " -n ", outputName, sep="")
        jobArgs <- append(jobArgs, systemCall)            
      }
    }
  }
  write(jobArgs, file=paste("RAXML_JobArgs", sep=""))
}
writejobArgs(files)

toRun <- scan(file="c70jobs.txt", what="character", sep="\n")
for(i in sequence(length(toRun))){
  system(toRun[i])
}




########################  

workingDir <- "~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLrepeatMLsearch"
setwd(workingDir)
library(ape)
library(phrynomics)

getBestTrees <- function(model, dataset){
  return(system(paste("ls RAxML_bestTree*", dataset, "*", model, "*", sep=""), intern=TRUE))
}

makeRFdistanceMatrix <- function(treeList){
  compare <- matrix(nrow=length(treeList), ncol=length(treeList))
  for(i in sequence(length(treeList))){
    tree1 <- read.tree(treeList[i])
    for(j in sequence(length(treeList))){
      if(i != j){
        tree2 <- read.tree(treeList[j])
        compare[i,j] <- phangorn::treedist(tree1, tree2)[[1]]  #symmetric difference
      }
    }
  }
return(compare)
}

models <- c("ASC", "GTR")
dataset <- paste("c", seq(from=5, to=70, by=5), "p3", sep="")
RFdists <- NULL
RFdistnames <- NULL
RFdistMatrix <- matrix(nrow=28, ncol=4)
place <- 0
for(i in sequence(length(models))){
  for(j in sequence(length(dataset))){
    place <- place+1
    print(place)
    RFdist <- NULL
    RFdistname <- NULL
    RFdist <- mean(makeRFdistanceMatrix(getBestTrees(models[i], dataset[j])), na.rm=TRUE)    
    RFdists <- c(RFdists, RFdist)
    RFdistname <- paste(models[i], dataset[j], sep="")
    print(paste(models[i], dataset[j], sep=""))
    RFdistnames <- c(RFdistnames, RFdistname)
    RFdistMatrix[place,] <- c(models[i], dataset[j], j, RFdist)
  }
 names(RFdists) <- RFdistnames
}

CheckInvocations <- function(workingDir){
  #read in info filenames - scrape out model, dataset, and rep from name
  #then grep for invocation line and scrape model and dataset to make sure it matches
  origWD <- getwd()
  setwd(workingDir)
  seeds <- NULL
  infoFiles <- system("ls *info*", intern=T)
  for(i in sequence(length(infoFiles))){
    modLine <- strsplit(infoFiles[i], "[.,_]")[[1]]
    if(any(grep("ASC", modLine)))
      mod1 <- "ASC_GTRCAT"
    else 
      mod1 <- "GTRCAT"
    dat1 <- strsplit(infoFiles[i], "[.,_]")[[1]][grep("\\d+", strsplit(infoFiles[i], "[.,_]")[[1]])[1]]
    dat1 <- gsub("noAmbigs", "", dat1)
    run <- strsplit(infoFiles[i], "[.,_]")[[1]][grep("\\d+", strsplit(infoFiles[i], "[.,_]")[[1]])[2]]
    invocLine <- strsplit(system(paste("grep -A 2 'RAxML was called as follows:'", infoFiles[i]), intern=TRUE)[3], split=" ")[[1]]
    mod2 <- invocLine[7]
    seeds <- c(seeds, invocLine[12])
    dat2 <- gsub("noAmbigs.snps", "", invocLine[5])
    if(mod1 != mod2)  
      print(paste("model invocation is wrong:", infoFiles[i]))
    if(dat1 != dat2)
      print(paste("dataset invocation is wrong:", infoFiles[i]))
  }
  setwd(origWD)
  return(as.numeric(seeds))
}

CheckSeeds <- function(seeds) {
  x1 <- NULL
  sorted.seeds <- sort(seeds)
  for (i in 1: length(seeds)-1) {
    x2 <- c(abs(sorted.seeds[i])-abs(sorted.seeds[i+1]))
    x1 <- c(x1, x2)
  }
  if ( 0 %in% x1)
    print(as.matrix(paste("repeat",  sort(seeds)[which(x1 == 0)])))
  if (!0 %in% x1)
    print("yahoo all are good!")
} 

systemCallSeeds <- CheckInvocations(workingDir)
CheckSeeds(systemCallSeeds)


# Print fig with two lines showing ASC mean RF vals and GTR mean RF vals
plot(as.numeric(RFdistMatrix[,3]), as.numeric(RFdistMatrix[,4]), xlab="dataset", ylab="Average RF", xaxt="n", type="n")
title(main="Average RF distances across 20 runs")
col1 <- "lightblue"
col2 <- "blue"
axis(side=1, at=1:14, labels=RFdistMatrix[1:14, 2])
ascsub <- RFdistMatrix[which(RFdistMatrix[,1] == "ASC"),]
gtrsub <- RFdistMatrix[which(RFdistMatrix[,1] == "GTR"),]
points(as.numeric(ascsub[,3]), as.numeric(ascsub[,4]), col=col1)
points(as.numeric(gtrsub[,3]), as.numeric(gtrsub[,4]), col=col2)
for(i in 1:13){
  segments(as.numeric(ascsub[i,3]), as.numeric(ascsub[i,4]), as.numeric(ascsub[i+1,3]), as.numeric(ascsub[i+1,4]), col=col1)
  segments(as.numeric(gtrsub[i,3]), as.numeric(gtrsub[i,4]), as.numeric(gtrsub[i+1,3]), as.numeric(gtrsub[i+1,4]), col=col2)
}
legtxt <- c("ASC", "GTR")
legcolors <- c(col1, col2)
legend("topleft", legend=legtxt, col=legcolors, lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Model"))) 




#Now to re analyze the TC/IC bootstrapped trees
#Here rerun and get new majority rule consensus tree for all runs
#then see if these improve support scatterplots

setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/TC-IC")

#raxmlHPC-PTHREADS-AVX -L MR -z <bootstrap trees> -m GTRCAT -n T1 -T4
BSfiles <- system("ls ../*bootstrap*", intern=TRUE)
rax <- "raxmlHPC-PTHREADS-SSE3"
for(i in sequence(length(BSfiles))){
  system(paste(rax, " -T 4 -L MR -z ", BSfiles[i], " -m GTRCAT -n ", gsub("RAxML_bootstrap.", "", gsub("../", "", BSfiles[i])), sep=""))
}
#scrape relTC and relTC-All off info files
infoFiles <- system("ls *info*", intern=TRUE)

GetTCIC <- function(infoFiles){
  TCICmat <- matrix(nrow=length(infoFiles), ncol=4)
  for(i in sequence(length(infoFiles))){
    splitName <- strsplit(infoFiles[i], "[.,_]")[[1]]
    TCICmat[i,1] <- splitName[4]
    TCICmat[i,2] <- splitName[3]
    lines <- system(paste("grep 'Relative tree certainty'", infoFiles[i]), intern=TRUE)
    TCICmat[i,3] <- strsplit(lines[1], "^\\D+", perl=TRUE)[[1]][2]
    TCICmat[i,4] <- strsplit(lines[2], "^\\D+", perl=TRUE)[[1]][2]
  }
colnames(TCICmat) <- c("model", "data", "TC", "TCA")
return(TCICmat)
}
TCIC <- GetTCIC(infoFiles)

orderToGo <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
plot(rep(1:14, 2), TCIC[,4], type="n", xlab="dataset", ylab="TC-ALL", xaxt="n")
axis(side=1, at=1:14, labels=orderToGo)
dataOrder <- rep(NA, 28)
for(i in 1:dim(TCIC)[1]){
  dataOrder[i] <- which(orderToGo == TCIC[i,2])
}
TCIC <- cbind(TCIC, dataOrder)
TCIC <- TCIC[order(as.numeric(TCIC[,5])),]
ascsub <- TCIC[which(TCIC[,1] == "ASC"),]
gtrsub <- TCIC[which(TCIC[,1] == "GTRCAT"),]
points(as.numeric(ascsub[,5]), as.numeric(ascsub[,4]), col=col1)
points(as.numeric(gtrsub[,5]), as.numeric(gtrsub[,4]), col=col2)
for(i in 1:13){
  segments(as.numeric(ascsub[i,5]), as.numeric(ascsub[i,4]), as.numeric(ascsub[i+1,5]), as.numeric(ascsub[i+1,4]), col=col1)
  segments(as.numeric(gtrsub[i,5]), as.numeric(gtrsub[i,4]), as.numeric(gtrsub[i+1,5]), as.numeric(gtrsub[i+1,4]), col=col2)
}
legtxt <- c("ASC", "GTR")
legcolors <- c(col1, col2)
legend("topright", legend=legtxt, col=legcolors, lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Model"))) 










































