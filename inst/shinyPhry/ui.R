##  --------------------------------  ##
##  ---      shinyPhrynomics     ---  ##
##  ---    written by B.Banbury  ---  ##
##  ---     v 1.1  27May2014   ---  ##
##  --------------------------------  ##


##  Load packages and source code
library(shiny)
library(ape)
library(phangorn)
library(devtools)
#devtools::install_github("bbanbury/phrynomics")
library(phrynomics)

##  shiny user interface
shinyUI(pageWithSidebar(  

##  ---          header          ---  ##
  headerPanel(
    HTML(paste("shiny", tags$span(style="color:Burlywood", "Phrynomics"), sep = ""))),



##  ---         side bar         ---  ##
  sidebarPanel(
    a("See phrynomics GitHub for full code base", href="https://github.com/bbanbury/phrynomics"),
    br(),
    a("exampleData.phy", href="example.txt"),
    br(),
    a("exampleData.nex", href="example.nex"),
    br(),
    h4("Original Dastset"),
    textInput("datasetName", "Enter SNP Dataset Name:", value="SNPdata"), 
    fileInput("SNPdataset", "Choose File To Upload (.phy, .nex, .txt, .snps):", accept=c(".snps", ".nex", ".txt", ".phy")),
    numericInput("obs", "Number Of Taxa to Preview:", 5),
    checkboxInput("rmInvSites", "Remove Invariant Sites", value=TRUE),
    checkboxInput("rmNonBin", "Remove Non-Binary Sites", value=TRUE),
    checkboxInput("takeRandom", "Take a single SNP from each locus (loci must be separated by a space)", value=FALSE),
    br(),
    h4("Transformations"),
    selectInput("transformChar", "Character for Missing Data:", choices=c("-", "?", "N"), selected="?"),
    h5("Transform to SNAPP"),
    h6(em("This will translate only binary sites. The most common allele will be translated to a 0, the least common to a 2, and the heterozygote to a 1.")),
    checkboxInput("transSNAPP", "Translate", value=FALSE),
    br(),
    h5("Transform to MrBayes Mk/Mkv Model"),
    h6(em("Will convert A,T,G,C to 1,2,3,4 (respectively), ambiguity codes are supported.")),
    checkboxInput("transMrBayes", "Translate", value=FALSE),
    br(),
    img(src="phryno.png", width=350),
    helpText("Photo by Jared Grummer"),
    verbatimTextOutput("phrynoversion")),

##  ---         main bar         ---  ##
  mainPanel(
    h4("Preview Original Dataset"),
    textOutput("OrigDataStats"),
    tableOutput("inputPreview"),
    h4("Preview Results"),
    textOutput("resultsDataStats"),
    tableOutput("resultsPreview"),
    verbatimTextOutput("communicationWindow"),
    hr(),
    downloadButton("downloadPhy", label="Download Full Results to Phy"),
    downloadButton("downloadNex", label="Download Full Results to Nexus"))
))
