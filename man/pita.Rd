\name{pita}
\alias{AddBLD}
\alias{AddTreeDist}
\alias{assTrees}
\alias{CalculateTotalTipBLError}
\alias{CheckSharedMonophy}
\alias{CompareMrBayesPosteriors}
\alias{CreateTreeList}
\alias{CreateTreeMatrix}
\alias{getColor}
\alias{GetCorresonding}
\alias{GetCorrespondingEdge}
\alias{GetEdgeList}
\alias{getMissingDataAmount}
\alias{GetMrBayesStatsPostAnalysis}
\alias{GetNumberOfSitesForLocus}
\alias{GetRAxMLStatsPostAnalysis}
\alias{GetSingleEdgeColor}
\alias{GetSpecies}
\alias{GetSupportValue}
\alias{MakeBranchLengthMatrix}
\alias{makeNumberwithcommas}
\alias{MakePresentAbsent}
\alias{MissingSpeciesVector}
\alias{nodeLeaves}
\alias{nodeOffspring}
\alias{ReturnMaxCI}
\alias{ReturnMinCI}

\title{A collection of pita tools for phrynomics that pertain to just us
}

\description{These probably won't be useful for anyone else
}

\usage{
AddBLD(TreeMatrixName, ListOfTrees)
AddTreeDist(TreeMatrixName, ListOfTrees)
assTrees(TreeMatrixName, ListOfTrees)
CalculateTotalTipBLError(BL.AllTrees)
CheckSharedMonophy(t1.node, tree1, tree2)
CompareMrBayesPosteriors(run1, run2, ntax=73)
CreateTreeList(filenames, analysis="RAxML")
CreateTreeMatrix(trees)
getColor(BLtable, nonSigColor="gray", sigColor="red", method=c("varOverlap", "meanWithin"))
GetCorresonding(corr.desc, t2)
GetCorrespondingEdge(t1.node, tree1, tree2)
GetEdgeList(tree)
getMissingDataAmount(filename)
GetMrBayesStatsPostAnalysis(workingDirectoryOfResults)
GetNumberOfSitesForLocus(data, locus)
GetRAxMLStatsPostAnalysis(workingDirectoryOfResults)
GetSingleEdgeColor(relativeBLdiff)
GetSpecies(taxa)
GetSupportValue(lineNumber, tstat.table)
MakeBranchLengthMatrix(tree1, tree2, analysis="RAxML", dataset=NULL)
makeNumberwithcommas(number)
MakePresentAbsent(dataset)
MissingSpeciesVector(data, SpeciesNames)
nodeLeaves(tree, anc.node)
nodeOffspring(tree, anc.node)
ReturnMaxCI(x, vstat)
ReturnMinCI(x, vstat)
}

\arguments{
  \item{TreeMatrixName}{Matrix of Tree names between ASC/GTR}
  \item{ListOfTrees}{List of Trees in the class phylo}
  \item{BL.AllTrees}{List of matrices of BLs}
  \item{t1.node}{Node of interest from tree 1}
  \item{tree1}{Tree 1}
  \item{tree2}{Tree 2}
  \item{run1}{Run 1 MrBayes--like ASC model}
  \item{run2}{Run 2 MrBayes--like GTR model}
  \item{ntax}{hard coded badness}
  \item{filenames}{vector of filenames of tree files}
  \item{analysis}{Which type of analysis, choices are "RAxML" or "MrBayes"}
  \item{trees}{all tree files}
  \item{BLtable}{Individual BL table}
  \item{nonSigColor}{which color represents the non-significant difference color}
  \item{sigColor}{which color represents the significant difference color}
  \item{method}{choices = c("varOverlap", "meanWithin")}
  \item{corr.desc}{descendent from tree1}
  \item{t2}{GetEdgeList matrix from tree2 info, this is the corresponding descendent}
  \item{t1.node}{GetEdgeList matrix from tree 1, node of interst}
  \item{tree}{A single phylogenetic tree}
  \item{filename}{For phrynomics named datasets (ex: "c70p3" will return 70)}
  \item{workingDirectoryOfResults}{path to working directory of MrBayes/RAxML runs}
  \item{data}{SNP dataset}
  \item{locus}{Which locus in the SNP dataset to calculate sites}
  \item{relativeBLdiff}{Relative BL difference}
  \item{taxa}{vector of taxa to extract species from (Ex: c("PH1", "PH2", "PH3") will return just "PH")}
  \item{linenumber}{some line number}
  \item{tstat.table}{some tstat table I don't remember}
  \item{number}{Any number}
  \item{SpeciesNames}{A vector of tip taxa}
  \item{anc.node}{ancestral node}
  \item{x}{don't remember}
  \item{vstat}{something having to do with vstat}  
}

\details{These are a collection of functions that we used for our paper, but probably need to be altered for the package.}
\value{returns a bunch of various stuffs...}
\author{BBanbury, but don't hold me to it}