setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/RAxML.autoBS/")

#scrape data off info files
#Get number bootsrtaps, get WRF ave splits, dataset

infoFiles <- system("ls *info*", intern=TRUE)
dataScrape <- data.frame(matrix(nrow=length(infoFiles), ncol=4))
colnames(dataScrape) <- c("model", "missing", "BS", "splits")
for(i in sequence(length(infoFiles))){
  dataScrape[i,1] <- strsplit(infoFiles[i], "[_,.]")[[1]][3]  #model
  dataScrape[i,2] <- strsplit(infoFiles[i], "\\D+")[[1]][2]  #missing data
  line <- system(paste("grep BS -A 1 ", infoFiles[i]), intern=TRUE)
  dataScrape[i,3] <- strsplit(line[1], "\\D+")[[1]][2]
  dataScrape[i,4] <- paste(strsplit(line[2], "\\D+")[[1]][3], strsplit(line[2], "\\D+")[[1]][4], sep=".")
}
dataScrape[,3] <- as.numeric(dataScrape[,3])
dataScrape[,2] <- as.numeric(dataScrape[,2])
dataScrape[,4] <- as.numeric(dataScrape[,4])
ASCdata <- dataScrape[which(dataScrape[,1] == "ASC"),]
GTRdata <- dataScrape[which(dataScrape[,1] == "GTR"),]

layout(matrix(1:2, nrow=1, byrow=TRUE), respect=TRUE)
orderToGo <- seq(from=5, to=70, by=5)
plot(dataScrape[,2], dataScrape[,3], type="n", xlab="missing data", ylab="number bootstraps")
for(i in 1:length(orderToGo)){
  points(ASCdata[which(ASCdata[,2] == orderToGo[i]),][,2], ASCdata[which(ASCdata[,2] == orderToGo[i]),][,3], col="blue")
  points(GTRdata[which(GTRdata[,2] == orderToGo[i]),][,2], GTRdata[which(GTRdata[,2] == orderToGo[i]),][,3], col="lightblue")
}
for(i in 1:(length(orderToGo)-1)){
  o1 <- orderToGo[i]
  o2 <- orderToGo[i+1]
  segments(ASCdata[which(ASCdata[,2] == o1),][,2], ASCdata[which(ASCdata[,2] == o1),][,3], ASCdata[which(ASCdata[,2] == o2),][,2], ASCdata[which(ASCdata[,2] == o2),][,3], col="blue")
  segments(GTRdata[which(GTRdata[,2] == o1),][,2], GTRdata[which(GTRdata[,2] == o1),][,3], GTRdata[which(GTRdata[,2] == o2),][,2], GTRdata[which(GTRdata[,2] == o2),][,3], col="lightblue")
}

plot(dataScrape[,2], dataScrape[,4], type="n", xlab="missing data", ylab="Ave Split")
for(i in 1:length(orderToGo)){
  points(ASCdata[which(ASCdata[,2] == orderToGo[i]),][,2], ASCdata[which(ASCdata[,2] == orderToGo[i]),][,4], col="blue")
  points(GTRdata[which(GTRdata[,2] == orderToGo[i]),][,2], GTRdata[which(GTRdata[,2] == orderToGo[i]),][,4], col="lightblue")
}
for(i in 1:(length(orderToGo)-1)){
  o1 <- orderToGo[i]
  o2 <- orderToGo[i+1]
  segments(ASCdata[which(ASCdata[,2] == o1),][,2], ASCdata[which(ASCdata[,2] == o1),][,4], ASCdata[which(ASCdata[,2] == o2),][,2], ASCdata[which(ASCdata[,2] == o2),][,4], col="blue")
  segments(GTRdata[which(GTRdata[,2] == o1),][,2], GTRdata[which(GTRdata[,2] == o1),][,4], GTRdata[which(GTRdata[,2] == o2),][,2], GTRdata[which(GTRdata[,2] == o2),][,4], col="lightblue")
}
