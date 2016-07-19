##  --------------------------------  ##
##  ---      shinyPhrynomics     ---  ##
##  ---    written by B.Banbury  ---  ##
##  ---     v 1.3  15August2014   ---  ##
##  --------------------------------  ##


##  Load packages and source code
library(shiny)
#devtools::install_github("shiny-incubator", "rstudio")

##  shiny user interface
shinyUI(pageWithSidebar(  

##  ---          header          ---  ##
  headerPanel(
    title=HTML(paste("shiny", tags$span(style="color:#BF00FF", "Phrynomics"), sep="")), windowTitle="shinyPhrynomics"),




##  ---         side bar         ---  ##
  sidebarPanel(
    a("See phrynomics GitHub for full code base", href="https://github.com/bbanbury/phrynomics"),
    br(),
    a("exampleData.phy", href="example.txt"),
    br(),
    checkboxInput("UseExampleData", "Click here to try it out with example data", value = FALSE),
    br(),
    br(),
    h4("Original Dastset"),
    textInput("datasetName", "Enter SNP Dataset Name:", value="SNPdata"), 
    fileInput("OrigData", "Choose File To Upload (.phy, .nex, .txt, .snps):", accept=c(".snps", ".nex", ".txt", ".phy")),
    textInput("minInds", "Enter Minimum number of individuals with sequence data required for a locus to be included in the dataset:", value=""), 
    br(),
    downloadButton("dRAxML", label="Transform and Download to RAxML"),
    downloadButton("dMrBayes", label="Transform and Download to MrBayes Mk/v Model"),
    downloadButton("dsnapp", label="Transform and Download to SNAPP"),
    br(), 
    br(), 
    br(), 
    img(src="phryno.png", width=350),
    helpText("Photo by Jared Grummer"),
    verbatimTextOutput("phrynoversion")),


##  ---         main bar         ---  ##
  mainPanel(
      h4("Vizualize"),
      br(),
      selectInput("plotIn", "Choose a Plot", choices=c("", "Per Locus", "Frequency", "Missing Data", "Heatmap"), selected=""),
      downloadButton("downloadSummPlot", label="Download Plot"),
      br(),
      img(plotOutput("summaryPlot1"))
  
    



#maybe make another tabPanel with information about us (NSF grant, paper, stats server, etc)
     
  ) #mainPanel

)) #shinyUI / pageWithSideBar
