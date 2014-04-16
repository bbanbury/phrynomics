###################################################
########## Load functions, trees, & data ##########
###################################################

library(ape)
library(phangorn)
source("~/phrynomics/trunk/phrynomicsFunctions.R")


# Load RAxML Trees (GTRGAMMA model) -- OLD
# analysis <- "RAxML"
# setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLResults")
# trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstrap support
# TreeList <- CreateTreeList(trees, "RAxML")
# treeMatrix <- CreateTreeMatrix(trees)

# Load RAxMLrun2 Trees (GTRCAT model)
analysis <- "RAxML"
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma")
trees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstrap support
TreeList <- CreateTreeList(trees, "RAxML")
treeMatrix <- CreateTreeMatrix(trees)

# Load MrBayes Trees
analysis <- "MrBayes"
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma")
trees <- system("ls *.con.tre", intern=T)  
TreeList <- CreateTreeList(trees, "MrBayes")
treeMatrix <- CreateTreeMatrix(trees)


setwd("~/Dropbox/UWstuff/phrynomics/Analyses/newFigs")


###################################################
######## Get whole-tree distance metrics ##########
###################################################

# Add whole-tree metrics (phangorn) data to the tree matrix
treeMatrix2 <- AddTreeDist(treeMatrix, TreeList)
# Add whole-tree metrics (Kuhner & Felsenstein; Kscore) data to tree matrix 
treeMatrix3 <- AddBLD(treeMatrix2, TreeList)


# Plot some whole tree measures
missingData <- as.numeric(sapply(rownames(treeMatrix3), function (x) strsplit(x, "\\D")[[1]][2]))
layout(matrix(1:3, nrow=1, byrow=TRUE), respect=TRUE)
plot(missingData, treeMatrix3$symmetric.difference, pch=21, bg="gray", main="Symmetric Difference")
plot(missingData, treeMatrix3$quadratic.path.difference, pch=21, bg="gray", main="Quadratic Path Difference")
plot(missingData, treeMatrix3$Kscore, pch=21, bg="gray", main="Kscore")



###################################################
## Compare Homologous branch lengths and support ##
###################################################

#setwd("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma")  #for MrBayes runs only
# Homologous Branch Comparisons between ASC & GTR trees
BL.AllTrees <- list()
for(i in sequence(dim(treeMatrix)[1])) {
  tree1 <- assTrees(treeMatrix3[i,1], TreeList)[[1]]
  tree2 <- assTrees(treeMatrix3[i,2], TreeList)[[1]]
  BL.AllTrees[[i]] <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset=rownames(treeMatrix)[i])
  names(BL.AllTrees)[[i]] <- rownames(treeMatrix)[i]
}
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/newFigs2")



# Plotting ASC trees with colored branches
orderToGo <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
RAxMLresults <- GetRAxMLStatsPostAnalysis("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma")
MrBayesresults <- GetMrBayesStatsPostAnalysis("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma")

for(i in sequence(length(orderToGo))){
  letters <- as.character(sequence(length(orderToGo)))
  #letters <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n")
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
  # nodelabels(c(NA, BL.AllTrees[[dataToUse]]$support[which(BL.AllTrees[[dataToUse]]$class == "internal")]), cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, -0.1))
  # nodelabels(c(NA, BL.AllTrees[[dataToUse]]$corr.support[which(BL.AllTrees[[dataToUse]]$class == "internal")]), cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, 1.1))
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



# Branch length density plots for paper
pdf(file=paste(analysis, "densityPlots.pdf", sep=""), width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(c(1:4), nrow=1, byrow=TRUE), respect=TRUE)
whichDatasets <- c("c5p3", "c25p3", "c45p3", "c65p3")
#whichDatasets <- c("c10p3", "c25p3", "c35p3", "c45p3")  #works for MrBayes trials too

for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == names(BL.AllTrees))
  BLs <- BL.AllTrees[[dataToUse]]$branchlength[which(BL.AllTrees[[dataToUse]]$present)]
  corr.BLs <- BL.AllTrees[[dataToUse]]$corr.BL[which(BL.AllTrees[[dataToUse]]$present)]
  maxy <- max(density(BLs)$y, density(corr.BLs)$y)
  maxx <- max(density(BLs)$x, density(corr.BLs)$x)
  plot(density(BLs), type="n", main="", frame=F, xlim=c(0, 0.35), ylim=c(0,160))
  lines(density(BLs), lwd=1, col="blue")
  lines(density(corr.BLs), lty=1, lwd=1)
  title(main=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""))
}
dev.off()

# Branch length density plots for supplement
pdf(file=paste(analysis, "all.densityPlots.pdf", sep=""), width=8.5, height=11)
op <- par(mar=par("mar")/1.7)
layout(matrix(c(1:16), nrow=4, byrow=TRUE), respect=TRUE)
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == names(BL.AllTrees))
  BLs <- BL.AllTrees[[dataToUse]]$branchlength[which(BL.AllTrees[[dataToUse]]$present)]
  corr.BLs <- BL.AllTrees[[dataToUse]]$corr.BL[which(BL.AllTrees[[dataToUse]]$present)]
  maxy <- max(density(BLs)$y, density(corr.BLs)$y)
  maxx <- max(density(BLs)$x, density(corr.BLs)$x)
  plot(density(BLs), type="n", main="", frame=F, xlim=c(0, 0.35), ylim=c(0,160))
  #leg.txt <- c("ASC", "GTR")
  #legend("topright", legend=leg.txt, text.col=c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)), merge = TRUE, bty="n", xjust=1, inset=0.02) 
  #legend("topright", legend=leg.txt, lwd=1, lty=1:2, merge = TRUE, bty="n", xjust=1, inset=0.02) 
  lines(density(BLs), lwd=1, col="blue")
  lines(density(corr.BLs), lty=1, lwd=1)
  #polygon(density(BLs), col=rgb(0,0,0,.3), border=NA)
  #polygon(density(corr.BLs), col=rgb(0,0,1,.3), border=NA)
  title(main=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""))
}
dev.off()

# Branch length scatter plots for paper
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
pdf(file=paste(analysis, "ScatterPlots.pdf", sep=""), width=8.5, height=5)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:4, nrow=1, byrow=TRUE), respect=TRUE)  
whichDatasets <- c("c5p3", "c25p3", "c45p3", "c65p3")
#whichDatasets <- c("c10p3", "c25p3", "c35p3", "c45p3")  #works for MrBayes trials too

for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == names(BL.AllTrees))
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
dev.off()

# Branch length scatter plots for supplement
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
pdf(file=paste(analysis, "all.ScatterPlots.pdf", sep=""), height=11, width=8.5)
op <- par(mar=par("mar")/1.7)
layout(matrix(c(1:16), nrow=4, byrow=TRUE), respect=TRUE)
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
#whichDatasets <- c("c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3")
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == names(BL.AllTrees))
  BLs <- BL.AllTrees[[dataToUse]]$branchlength[which(BL.AllTrees[[dataToUse]]$present)]
  corr.BLs <- BL.AllTrees[[dataToUse]]$corr.BL[which(BL.AllTrees[[dataToUse]]$present)]
  colors <- getColor(BL.AllTrees[[dataToUse]], method="meanWithin")[which(BL.AllTrees[[dataToUse]]$present)]
  plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="non-ASC", type="n")
  linmod <- lm(corr.BLs ~ BLs)
  lines(c(-1,1), c(-1,1))
  abline(linmod, lty=2)
  y <- 0.18
  points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
  points(BLs, corr.BLs, pch=21, bg=colors)
  text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
  title(main=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""))
}
dev.off()

# Make Scatter plots for support
#MrBayes needs to be in the dir: setwd("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma")
pdf(file=paste(analysis, "all.ScatterPlots.support.pdf", sep=""), width=8.5, height=11)
op <- par(mar=par("mar")/1.7)
layout(matrix(1:16, nrow=4, byrow=TRUE), respect=TRUE)
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
if(analysis == "RAxML"){
  for(i in sequence(dim(treeMatrix)[1])){
    dataToUse <- which(whichDatasets[i] == names(BL.AllTrees))
    datarows <- which(BL.AllTrees[[dataToUse]]$present)[which(BL.AllTrees[[dataToUse]]$present) %in% which(BL.AllTrees[[dataToUse]]$support != 0)]  #don't want all the tip support 0s or non homologous clades
    support <- BL.AllTrees[[dataToUse]]$support[datarows]
    corr.support <- BL.AllTrees[[dataToUse]]$corr.support[datarows]
    if(analysis == "RAxML")
      xlims <- ylims <- c(1,100)
    if(analysis == "MrBayes")
      xlims <- ylims <- c(0,1)
    plot(support, corr.support, ylab="GTR Tree", xlab="ASC Tree", pch=21, bg="gray", xlim=xlims, ylim=ylims)
    title(main=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""))
  }
  dev.off()
}
if(analysis == "MrBayes"){
  for(i in sequence(length(whichDatasets))){
    dataToUse <- which(whichDatasets[i] == rownames(treeMatrix))
    supportMatrix <- CompareMrBayesPosteriors(treeMatrix[dataToUse,1], treeMatrix[dataToUse,2])
    xlims <- ylims <- c(0,1)
    plot(supportMatrix$ASC.support, supportMatrix$GTR.support, ylab="GTR Tree", xlab="ASC Tree", pch=21, bg="gray", xlim=xlims, ylim=ylims)
    title(main=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""))
  }
  dev.off()
}



# Get Taxa SNP overlap 
datasets <- system("ls ~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma/c*.snps", intern=T)
dataOverlap <- list()
for(i in sequence(length(datasets))){
  dataset <- read.table(datasets[[i]], row.names=1, colClasses="character", skip=1)
  dataOverlap[[i]] <- DataOverlap(dataset)[[2]]
  names(dataOverlap)[[i]] <- strsplit(strsplit(datasets[i], "/")[[1]][10], "no")[[1]][1]
}

# Get Taxa total branchlengh differences from root to tip
sum.BLdiff <- list()
for(i in sequence(length(BL.AllTrees))){
  sum.BLdiff[[i]] <- CalculateTotalTipBLError(BL.AllTrees[[i]])
  names(sum.BLdiff)[[i]] <- names(BL.AllTrees)[[i]]
}

#Get Taxon specific branchlength error
taxon.BLdiff <- list()
for(i in sequence(length(BL.AllTrees))){
  taxon.BLdiff[[i]] <- GetJustTipBLError(BL.AllTrees[[i]])
  names(taxon.BLdiff)[[i]] <- names(BL.AllTrees)[[i]]
}

#plot just data overlap
layout(matrix(c(1:16), nrow=4, byrow=TRUE), respect=TRUE)
for(i in sequence(length(dataOverlap))){
  whichData <- names(dataOverlap[i])
  plot(dataOverlap[[which(names(dataOverlap) == whichData)]], ylab="data Overlap")
  title(sub=whichData)
}

# sum of BL differences
layout(matrix(c(1:16), nrow=4, byrow=TRUE), respect=TRUE)
for(i in sequence(length(dataOverlap))){
  orderToGo <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
  whichData <- orderToGo[i]
  plot(dataOverlap[[which(names(dataOverlap) == whichData)]], sum.BLdiff[[which(names(sum.BLdiff) == whichData)]], xlab="data Overlap", ylab="sum BL differences")
  title(sub=whichData)
}

# data overlap with tip BL differences
layout(matrix(c(1:16), nrow=4, byrow=TRUE), respect=TRUE)
for(i in sequence(length(dataOverlap))){
 # orderToGo <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
  orderToGo <- c("c5p3", "c25p3", "c45p3", "c65p3")

  whichData <- orderToGo[i]
  plot(dataOverlap[[which(names(dataOverlap) == whichData)]], taxon.BLdiff[[which(names(taxon.BLdiff) == whichData)]], xlab="data Overlap", ylab="taxon BL Diff", type="n")
  whichPhrynosoma <- grep(pattern="PH", names(dataOverlap[[which(names(dataOverlap) == whichData)]]))
  points(dataOverlap[[which(names(dataOverlap) == whichData)]][whichPhrynosoma], taxon.BLdiff[[which(names(taxon.BLdiff) == whichData)]][whichPhrynosoma], col="red")
  points(dataOverlap[[which(names(dataOverlap) == whichData)]][-whichPhrynosoma], taxon.BLdiff[[which(names(taxon.BLdiff) == whichData)]][-whichPhrynosoma], col="blue")
  title(sub=whichData)
}

#boxplot the data
#make a new table with all data
all.data <- NULL
orderToGo <- c("c5p3", "c25p3", "c45p3", "c65p3")
for(i in sequence(length(orderToGo))){
  whichData <- orderToGo[i]
  data1 <- taxon.BLdiff[[which(names(taxon.BLdiff) == whichData)]]
  data2 <- rep("OG", length(data1))
  data2[grep(pattern="PH", names(dataOverlap[[which(names(dataOverlap) == whichData)]]))] <- "PH"
  data2 <- paste(whichData, data2, sep=".")
  data3 <- data.frame(data1, data2)
  colnames(data3) <- c("BL", "whichData") 
  all.data <- rbind(all.data, data3)
}
boxplot(all.data[,1] ~ all.data[,2], range=0)









###################################################
######## Scrape Data from post-analysis ###########
###################################################

ML.results <- GetRAxMLStatsPostAnalysis("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma")
MrB.results <- MrBayesresults

#results <- ML.results

# Plot number alignment patterns vs missing data
plot(ML.results$AlignmentPatterns[which(ML.results$Model == "ASC")], ML.results$MissingData[which(ML.results$Model == "ASC")], type="n", ylab="Percent Missing Data", xlab="Distinct Alignment Patterns")
text(ML.results$AlignmentPatterns[which(ML.results$Model == "ASC")], ML.results$MissingData[which(ML.results$Model == "ASC")], labels=ML.results$Level[which(ML.results$Model == "ASC")])

# Plot likelihood scores by model
plot(ML.results$Likelihood[which(ML.results$Model == "ASC")], ML.results$Likelihood[which(ML.results$Model == "GTRCAT")], xlab="ASC Likelihood", ylab="GTR Likelihood")
linMod <- lm(ML.results$Likelihood[which(ML.results$Model == "ASC")] ~ ML.results$Likelihood[which(ML.results$Model == "GTR")], ML.results)
abline(linMod)
cor(ML.results$Likelihood[which(ML.results$Model == "ASC")], ML.results$Likelihood[which(ML.results$Model == "GTRCAT")])
text(-100000, 100, paste("r =", cor(ML.results$Likelihood[which(ML.results$Model == "ASC")], ML.results$Likelihood[which(ML.results$Model == "GTRCAT")])), pos=2) 

# Plot bootstrap times by model and missing data
levels <- seq(from=5, to=70, by=5)
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
plot(levels, ML.results$BootstrapTime[which(ML.results$Model == "ASC")], type="n")
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  points(levels[i], ML.results$BootstrapTime[dataToUse[1]], pch=21, bg="lightblue")
  points(levels[i], ML.results$BootstrapTime[dataToUse[2]], pch=21, bg="blue")
}
legend("topright", legend=c("GTR", "GTR+ASC"), col=c("blue", "lightblue"), lwd=1, merge = TRUE, bty="n", xjust=1, inset=0.05, cex=1) 



# Plot alpha by levels
levels <- seq(from=5, to=70, by=5)
MrB.levels <- seq(from=10, to=45, by=5)
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
MrBdatasets <-  c("c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3")
layout(matrix(1:2, nrow=1, byrow=TRUE), respect=TRUE)
plot(c(levels, levels), ML.results$Alpha, type="n", xlab="Missing Data", ylab="Alpha (Gamma Shape Parameter)")
title(main="RAxML")
legend("bottomleft", legend=c("GTR", "GTR+ASC"), col=c("blue", "lightblue"), lwd=1, merge = TRUE, bty="n", xjust=1, inset=0.05, cex=1) 
for(i in sequence(length(whichDatasets)-1)){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  nextDataToUse <- which(whichDatasets[i+1] == ML.results$Level)
  segments(levels[i], ML.results$Alpha[dataToUse[1]], levels[i+1], ML.results$Alpha[nextDataToUse[1]], col="lightblue")
  segments(levels[i], ML.results$Alpha[dataToUse[2]], levels[i+1], ML.results$Alpha[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  points(levels[i], ML.results$Alpha[dataToUse[1]], pch=21, bg="lightblue")
  points(levels[i], ML.results$Alpha[dataToUse[2]], pch=21, bg="blue")
}
plot(c(MrB.levels, MrB.levels), MrB.results$Alpha, type="n", xlab="Missing Data", ylab="Alpha (Gamma Shape Parameter)")
title(main="MrBayes")
legend("bottomleft", legend=c("GTR", "GTR+ASC"), col=c("blue", "lightblue"), lwd=1, merge = TRUE, bty="n", xjust=1, inset=0.05, cex=1) 
for(i in sequence(length(MrBdatasets)-1)){
  dataToUse <- which(MrBdatasets[i] == MrB.results$Level)
  nextDataToUse <- which(MrBdatasets[i+1] == MrB.results$Level)
  segments(MrB.levels[i], MrB.results$Alpha[dataToUse[1]], MrB.levels[i+1], MrB.results$Alpha[nextDataToUse[1]], col="lightblue")
  segments(MrB.levels[i], MrB.results$Alpha[dataToUse[2]], MrB.levels[i+1], MrB.results$Alpha[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(MrBdatasets))){
  dataToUse <- which(MrBdatasets[i] == MrB.results$Level)
  if(!is.null(MrB.results$Alpha.lowCI)){
    arrows(MrB.levels[i], MrB.results$Alpha.lowCI[dataToUse[1]], MrB.levels[i], MrB.results$Alpha.uppCI[dataToUse[1]], code=3, length=0.05, col="lightblue", angle=90)
    arrows(MrB.levels[i], MrB.results$Alpha.lowCI[dataToUse[2]], MrB.levels[i], MrB.results$Alpha.uppCI[dataToUse[2]], code=3, length=0.05, col="blue", angle=90)
  }
  points(MrB.levels[i], MrB.results$Alpha[dataToUse[1]], pch=21, bg="lightblue")
  points(MrB.levels[i], MrB.results$Alpha[dataToUse[2]], pch=21, bg="blue")
}



#Plot tree length by amount of missing data
levels <- seq(from=5, to=70, by=5)
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
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
plot(c(levels, levels), ML.results$TreeLength, type="n", ylab="Tree Length", xlab="Missing Data")
title(main="MrBayes")
legend("topright", legend=c("GTR", "GTR + ASC"), col=c("blue", "lightblue"), lwd=1, merge=TRUE, bty="n", xjust=1, inset=0.02, cex=1) 
for(i in sequence(length(whichDatasets)-1)){
  dataToUse <- which(whichDatasets[i] == MrB.results$Level)
  nextDataToUse <- which(whichDatasets[i+1] == MrB.results$Level)
  segments(levels[i], MrB.results$TreeLength[dataToUse[1]], levels[i+1], MrB.results$TreeLength[nextDataToUse[1]], col="lightblue")
  segments(levels[i], MrB.results$TreeLength[dataToUse[2]], levels[i+1], MrB.results$TreeLength[nextDataToUse[2]], col="blue")  
}
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == MrB.results$Level)
  arrows(levels[i], MrB.results$TreeLength.lowCI[dataToUse[1]], levels[i], MrB.results$TreeLength.uppCI[dataToUse[1]], code=3, length=0.05, col="lightblue", angle=90)
  arrows(levels[i], MrB.results$TreeLength.lowCI[dataToUse[2]], levels[i], MrB.results$TreeLength.uppCI[dataToUse[2]], code=3, length=0.05, col="blue", angle=90)
  points(levels[i], MrB.results$TreeLength[dataToUse[1]], pch=21, bg="lightblue")
  points(levels[i], MrB.results$TreeLength[dataToUse[2]], pch=21, bg="blue")
}




# Compare tree lengths by model
whichDatasets <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
plot(ML.results$TreeLength[which(ML.results$Model == "ASC")], ML.results$TreeLength[which(ML.results$Model == "GTRCAT")], type="n", xlab="ASC Tree Length", ylab="GTR Tree Length", xlim=c(0,7))
for(i in sequence(length(whichDatasets)-1)){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  nextDataToUse <- which(whichDatasets[i+1] == ML.results$Level)
  segments(ML.results$TreeLength[dataToUse[1]], ML.results$TreeLength[dataToUse[2]], ML.results$TreeLength[nextDataToUse[1]], ML.results$TreeLength[nextDataToUse[2]], col="blue")
}
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == ML.results$Level)
  points(ML.results$TreeLength[dataToUse[1]], ML.results$TreeLength[dataToUse[2]], col="blue", cex=3, bg="white", pch=21)  
  text(ML.results$TreeLength[dataToUse[1]], ML.results$TreeLength[dataToUse[2]], labels=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""), col="blue", cex=0.75)
}
for(i in sequence(length(whichDatasets)-1)){
  dataToUse <- which(whichDatasets[i] == MrB.results$Level)
  nextDataToUse <- which(whichDatasets[i+1] == MrB.results$Level)
  segments(MrB.results$TreeLength[dataToUse[1]], MrB.results$TreeLength[dataToUse[2]], MrB.results$TreeLength[nextDataToUse[1]], MrB.results$TreeLength[nextDataToUse[2]], col="green")
}
for(i in sequence(length(whichDatasets))){
  dataToUse <- which(whichDatasets[i] == MrB.results$Level)
  points(MrB.results$TreeLength[dataToUse[1]], MrB.results$TreeLength[dataToUse[2]], col="green", cex=3, bg="white", pch=21)  
  text(MrB.results$TreeLength[dataToUse[1]], MrB.results$TreeLength[dataToUse[2]], labels=paste("s", strsplit(whichDatasets[[i]], "\\D")[[1]][2], sep=""), col="green", cex=0.75)
}




############################################
###  Examine the effect of using all  ######
###  sites, not just variables   ###########
############################################

#compare branchlengths for ASC_c55p3 and GTR_allSites
ASCtree <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns2/RAxML_bipartitions.c55p3_ASC_GTRCAT")
allSitestree <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.c55p3_allSites")
tree1 <- ASCtree
tree2 <- allSitestree
tree2$tip.label[33] -> tree1$tip.label[46]  #different tip labels
BL.AllTrees <- MakeBranchLengthMatrix(tree1, tree2, analysis=analysis, dataset="c55p3")

#branchlengths
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
pdf(file="BL.ASC_vs_allSites", width=5, height=5)
BLs <- BL.AllTrees$branchlength[which(BL.AllTrees$present)]
corr.BLs <- BL.AllTrees$corr.BL[which(BL.AllTrees$present)]
plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="non-ASC", type="n")
linmod <- lm(corr.BLs ~ BLs)
abline(linmod, lty=2)
y <- 0.18
points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
points(BLs, corr.BLs, pch=21, bg="gray")
text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
lines(c(-1,1), c(-1,1))
title(main=paste("s55"))
dev.off()

#density branchlength
pdf(file="densityPlots.pdf", width=8.5, height=5)
BLs <- BL.AllTrees$branchlength[which(BL.AllTrees$present)]
corr.BLs <- BL.AllTrees$corr.BL[which(BL.AllTrees$present)]
maxy <- max(density(BLs)$y, density(corr.BLs)$y)
maxx <- max(density(BLs)$x, density(corr.BLs)$x)
plot(density(BLs), type="n", main="", frame=F, xlim=c(0, 0.35), ylim=c(0,160))
lines(density(BLs), lwd=1, col="blue")
lines(density(corr.BLs), lty=1, lwd=1)
dev.off()



#color tree
pdf(file="ASCvsallSites.pdf", width=8.5, height=11)
plot(tree1, edge.lty=BL.AllTrees$edgelty, edge.color=BL.AllTrees$edgeColor, cex=0.5, edge.width=2)
legtxt <- c("Discordant", "< -10%", "-10% to 10%", "> 10%", "> 20%", "> 30%", "> 40%", "> 50%")
legcolors <- c("gray", rgb(51,51,255, max=255), "gray", rgb(255,255,102, max=255), rgb(255,178,102, max=255), rgb(225,128,0, max=255), rgb(225,0,0, max=255), rgb(153,0,0, max=255))
legend("bottomleft", legend=legtxt, col=legcolors, lty=c(2,rep(1,7)), lwd=2, merge = TRUE, bty="n", xjust=1, inset=0.02, cex=0.75, title=expression(underline("Branchlength Difference"))) 
nodelabels(c(NA, BL.AllTrees$support[which(BL.AllTrees$class == "internal")]), cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, -0.1))
nodelabels(c(NA, BL.AllTrees$corr.support[which(BL.AllTrees$class == "internal")]), cex=0.5, col="black", bg="white", frame="none", adj=c(1.1, 1.1))
dev.off()


#compare branchlengths of three c55 trees (full, nonASC, and ASC) with made up invariant sites
full <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.full.c55")  #dataset=c55.extraData
nonASC <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.nonASC.c55")  #dataset=c55.original
ASC <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.ASC.c55")   #dataset=c55.original
full2 <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.full2.c55")
full3 <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.full3.c55")
BL.AllTrees1 <- MakeBranchLengthMatrix(ASC, nonASC, analysis="RAxML", dataset="c55p3")
BL.AllTrees2 <- MakeBranchLengthMatrix(ASC, full, analysis="RAxML", dataset="c55p3")
BL.AllTrees3 <- MakeBranchLengthMatrix(ASC, full2, analysis="RAxML", dataset="c55p3")
BL.AllTrees4 <- MakeBranchLengthMatrix(ASC, full3, analysis="RAxML", dataset="c55p3")
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x

layout(matrix(1:4, nrow=1, byrow=TRUE), respect=TRUE)
BLs <- BL.AllTrees1$branchlength[which(BL.AllTrees1$present)]
corr.BLs <- BL.AllTrees1$corr.BL[which(BL.AllTrees1$present)]
plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="nonASC", type="n")
linmod <- lm(corr.BLs ~ BLs)
abline(linmod, lty=2)
y <- 0.18
points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
points(BLs, corr.BLs, pch=21, bg="gray")
text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
lines(c(-1,1), c(-1,1))
title(main=paste("s55"))

BLs <- BL.AllTrees2$branchlength[which(BL.AllTrees2$present)]
corr.BLs <- BL.AllTrees2$corr.BL[which(BL.AllTrees2$present)]
plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="full", type="n")
linmod <- lm(corr.BLs ~ BLs)
abline(linmod, lty=2)
y <- 0.18
points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
points(BLs, corr.BLs, pch=21, bg="gray")
text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
lines(c(-1,1), c(-1,1))
title(main=paste("s55"))

BLs <- BL.AllTrees3$branchlength[which(BL.AllTrees3$present)]
corr.BLs <- BL.AllTrees3$corr.BL[which(BL.AllTrees3$present)]
plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="full2", type="n")
linmod <- lm(corr.BLs ~ BLs)
abline(linmod, lty=2)
y <- 0.18
points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
points(BLs, corr.BLs, pch=21, bg="gray")
text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
lines(c(-1,1), c(-1,1))
title(main=paste("s55"))

BLs <- BL.AllTrees4$branchlength[which(BL.AllTrees4$present)]
corr.BLs <- BL.AllTrees4$corr.BL[which(BL.AllTrees4$present)]
plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, 0.22), xlim=c(0, 0.22), xlab="ASC", ylab="full+gamma", type="n")
linmod <- lm(corr.BLs ~ BLs)
abline(linmod, lty=2)
y <- 0.18
points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
points(BLs, corr.BLs, pch=21, bg="gray")
text(x=getX(y, linmod), y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
lines(c(-1,1), c(-1,1))
title(main=paste("s55"))


#For simulated data, check branchlengths of 1) ASC v nonASC, 2) ASC v full, 3) ASC v ifull
full <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.c45.full")
gfull <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.c45.gfull")
igfull <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.c45.igfull")
nonASC <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.c45.nonASC")
ASC <- read.tree("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites/RAxML_bipartitions.c45.ASC")

BL.AllTrees <- list()
BL.AllTrees[[1]] <- MakeBranchLengthMatrix(ASC, nonASC, analysis="RAxML", dataset="c45")
BL.AllTrees[[2]] <- MakeBranchLengthMatrix(ASC, full, analysis="RAxML", dataset="c45")
BL.AllTrees[[3]] <- MakeBranchLengthMatrix(ASC, gfull, analysis="RAxML", dataset="c45")
BL.AllTrees[[4]] <- MakeBranchLengthMatrix(ASC, igfull, analysis="RAxML", dataset="c45")
names(BL.AllTrees) <- c("nonASC", "full", "gfull", "igfull")

#make BL scatter plots
usr2dev <- function(x) 180/pi*atan(x*diff(par("usr")[1:2])/diff(par("usr")[3:4]))
getX <- function(y, linmod) (y-linmod$coefficients[1])/linmod$coefficients[2]  #(y -b)/m = x
layout(matrix(1:length(BL.AllTrees), nrow=1, byrow=TRUE), respect=TRUE)
for(i in sequence(length(BL.AllTrees))){
  BLs <- BL.AllTrees[[i]]$branchlength[which(BL.AllTrees[[i]]$present)]
  corr.BLs <- BL.AllTrees[[i]]$corr.BL[which(BL.AllTrees[[i]]$present)]
  maxBL <- max(c(BLs, corr.BLs))
  plot(BLs, corr.BLs, pch=21, bg="gray", ylim=c(0, maxBL), xlim=c(0, maxBL), xlab="ASC", ylab=names(BL.AllTrees)[[i]], type="n")
  linmod <- lm(corr.BLs ~ BLs)
  abline(linmod, lty=2)
  y <- maxBL - (maxBL*.33)
  x <- getX(y, linmod) - (getX(y, linmod) * .33)
  points(getX(y, linmod), y, col="white", pch=21, bg="white", cex=12, crt=usr2dev(linmod$coefficients[2]))
  points(BLs, corr.BLs, pch=21, bg="gray")
  text(x=x, y=y, paste("m =", round(linmod$coefficients[2], digits=2)), srt=usr2dev(linmod$coefficients[2]))
  lines(c(-1,1), c(-1,1))
  title(main=paste("ASC v", names(BL.AllTrees)[[i]]))
}







############################################
###  Get a list of trees for TreeViz  ######
############################################

# Load RAxMLrun2 Trees (GTRCAT model)
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma")
RAxMLtrees <- system("ls RAxML_bipartitions.*", intern=T)  #trees with bootstrap support
TreeList1 <- CreateTreeList(RAxMLtrees, "RAxML")
treeMatrix1 <- CreateTreeMatrix(RAxMLtrees)

# Load MrBayes Trees
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma")
MrBayestrees <- system("ls *.con.tre", intern=T)  
TreeList2 <- CreateTreeList(MrBayestrees, "MrBayes")
treeMatrix2 <- CreateTreeMatrix(MrBayestrees)

setwd("~/Dropbox/UWstuff/phrynomics/Analyses")
write("", file="NewickTrees.txt")
orderToGo <- c("c5p3", "c10p3", "c15p3", "c20p3", "c25p3", "c30p3", "c35p3", "c40p3", "c45p3", "c50p3", "c55p3", "c60p3", "c65p3", "c70p3")
ns <- 0
for(m in 1:dim(treeMatrix1)[2]){
  for(i in sequence(length(orderToGo))) {
    treeToUse <- treeMatrix1[which(rownames(treeMatrix1) == orderToGo[i]), m]
    print(treeToUse)
    tree <- TreeList1[which(names(TreeList1) == treeToUse)][[1]]
    tree$node.label <- NULL
    write.tree(tree, file="NewickTrees.txt", append=TRUE)
    ns <- ns+1
  }
}
for(m in 1:dim(treeMatrix2)[2]){
  for(i in sequence(length(orderToGo))) {
    treeToUse <- treeMatrix2[which(rownames(treeMatrix2) == orderToGo[i]), m]
    print(treeToUse)
    tree <- TreeList2[which(names(TreeList2) == treeToUse)][[1]]
    tree$node.label <- NULL
    write.tree(tree, file="NewickTrees.txt", append=TRUE)
    ns <- ns+1
  }
}







RAxMLresults <- GetRAxMLStatsPostAnalysis("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxMLruns.noGamma")
MrBayesresults <- GetMrBayesStatsPostAnalysis("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns/MrBayes.noGamma")















