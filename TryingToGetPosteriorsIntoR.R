file <- "tree2.nex"
#file <- "~/Dropbox/Barb/5GBAYES.nex.con.tre"
#file <- "~/testMrBayes.nex.con.tre"

GetTreeLine <- function(file){
#for single trees ONLY
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  X <- X[grep("tree", X, ignore.case=TRUE)][-1]  #grep fro "tree" and then take second pices, which should be the tree (first is "begin trees" block)
  if(length(X) > 1)
    stop("There is more than one tree here and we can only do one for now")
  X <- paste(strsplit(X, "")[[1]][-1:-gregexpr("]", X)[[1]][1]], collapse="")
  X <- gsub("^\\s+", "", X)
  return(X)
}
#GetTreeLine(file)

GetStats <- function (file) {
  initTree <- read.nexus(file)
  TreeLength <- length(initTree$tip.label) + initTree$Nnode
  X <- GetTreeLine(file)
  test <- gsub("\\[[^]]*\\]", "", X)
  tab <- unlist(strsplit(X, "\\["))[-1] 
  if(length(tab) != TreeLength)
    stop("problem here") 
  tab <- gsub("&|;|\\]", "", tab)
  tab <- gsub(":.+$", "", tab)
  foo <- function(x) {
    return(unlist(strsplit(x, ",")))
  }
  tab <- lapply(tab, foo)
  for (i in seq(along = tab)) {
    ind <- grep("[{]", tab[[i]])
    names <- gsub("=.+$", "", tab[[i]][ind])
    tab[[i]][ind] <- gsub("[{]", "", tab[[i]][ind])
    tab[[i]][ind] <- gsub("=", "_MIN=", tab[[i]][ind])
    tab[[i]][ind + 1] <- gsub("[}]", "", tab[[i]][ind + 1])
    tab[[i]][ind + 1] <- paste(paste(names, "MAX=", sep = "_"), tab[[i]][ind + 1])
  }
  ttab <- data.frame()
  stats <- unique(gsub("=.+$", "", unlist(tab)))
  for (i in seq(along = tab)) {
    for (j in seq(along = stats)) {
      ind <- grep(paste(stats[j], "=", sep = ""), tab[[i]], fixed=T)
      if (length(ind) > 0) {
        v <- strsplit(gsub(paste(stats[j], "=", sep = ""), "", tab[[i]][ind], fixed=T), "")[[1]]
        whichNumbers <- grep("\\d+", v)
        v <- paste(v[min(whichNumbers):max(whichNumbers)], collapse="")
        ttab[i, j] <- v
      }
    }
  }
  colnames(ttab) <- stats
  tip <- which(is.na(ttab$prob))
  return(ttab)
}
GetStats(file)


read.figtree <- function (file, nodeLabelsToKeep="prob", digits = NULL) {
  #if(!is.binary.tree(read.nexus(file)))
  #  stop("Tree is not binary so you can't do this yet...")
  tab <- GetStats(file)
  if(!nodeLabelsToKeep %in% colnames(tab))
    stop(paste("nodeLabelsToKeep must be one of", paste(colnames(tab), collapse=", ")))
  X <- GetTreeLine(file)
  interior <- which(!is.na(tab$prob))
  #interior <- c(3,7,8,9)  
  tree <- gsub("\\[[^]]*\\]", "", X)
  brl <- unlist(strsplit(tree, ":"))[-1]
  brl <- gsub("[( | ) | ;]", "", brl)
  brl <- strsplit(brl, ",")
  foo <- function(x) x <- head(x, 1)
  brl <- unlist(lapply(brl, foo))
  brl <- paste("", brl, sep = ":")
  nodestats <- vector(mode = "list", length = dim(tab)[2])
  for (i in seq(along = nodestats)) {
    newtree <- tree
    newtree <- gsub(";", "", newtree)
    val <- tab[, i]
    ggg <- paste(val, brl, sep = "")
    ggg[length(ggg)] <- paste(tail(val, 1), ";", sep = "")
    a <- strsplit(newtree, ":")[[1]]
    for(j in sequence(length(interior))){  
      a[interior[j]] <- paste(a[interior[j]], val[interior[j]], sep="")
    }
    newtree <- paste(paste(a, collapse=":"), ";", sep="")
    dt <- read.tree(text = newtree)
    z <- dt$node.label
    nodestats[[i]] <- z
    names(nodestats)[i] <- colnames(tab)[i]
  }
  tr <- read.nexus(file)
  whichStat <- which(names(nodestats) == nodeLabelsToKeep)
  tr$node.label <- nodestats[whichStat][[1]]
  class(tr) <- ("phylo")
  attr(tr, "origin") <- file
  return(tr)
}
read.figtree(file)






