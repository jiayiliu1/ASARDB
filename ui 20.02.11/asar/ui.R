#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(d3heatmap)
library(RColorBrewer)
library(rhandsontable)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  titlePanel(maintitle),
  sidebarPanel(
    selectInput("normMethod", "Select Normalization Method:",
                choices = c("None" = "none", "TMM" = "TMM", "DESeq2" = "DESeq2", "RPKM" = "RPKM")),
    conditionalPanel(condition = "input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     uiOutput("Mgall")),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     uiOutput("mg1")),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4 || input.conditionedPanels==5",
                     selectInput(inputId = "tl1", label = taxone, choices = tax1n, selected = tax1selected, selectize = FALSE),
                     uiOutput("taxNames")),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==3",
                     selectInput(inputId = "tl2", label = taxtwo, choices = tax2n, selected = tax2selected)),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2 || input.conditionedPanels==3",
                     selectInput(inputId = "fl1", label = funcone, choices = func1n, selected = func1selected, selectize = FALSE),
                     uiOutput("funNames")),
    conditionalPanel(condition = "input.conditionedPanels==1 || input.conditionedPanels==2",
                     selectInput(inputId = "fl2", label = functwo, choices = func2n, selected = func2selected)),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     sliderInput("pix1", "Heatmap height (px)", value = 20, min = 10, max = 100)),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     sliderInput("pix2", "Heatmap height (px)", value = 20, min = 10, max = 100)),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     sliderInput("pix3", "Heatmap height (px)", value = 20, min = 10, max = 100)),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     actionButton("goButton", "GO"),
                     sliderInput("pix4", "Heatmap height (px)", value = 20, min = 10, max = 100)
    ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("path", "GO"),
                     uiOutput("PathwayID")
    ),
################ Truncation conditional pannel ############################
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3",
                     #selectInput(inputId = "prfun", label = pruneone, choices = prunec1n, selected = prunec1selected, selectize = FALSE),
                     radioButtons(inputId = "trfun", label = truncone, choices = list("sum", "max","sd"),inline=TRUE),
                     textInput("trnum","Rows:",value = 50)
    ),
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
                     textInput("filename","Enter file name")
    ),
    conditionalPanel(condition = "input.conditionedPanels==1 ||input.conditionedPanels==2 || input.conditionedPanels==3 || input.conditionedPanels==4",
                     radioButtons(inputId = "var3", label = "Select the file format", choices = list("png", "pdf"))
    ),
    conditionalPanel(condition = "input.conditionedPanels==1",
                     uiOutput("downLink1")
    ),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     uiOutput("downLink2")
    ),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     uiOutput("downLink3")
    ),
    conditionalPanel(condition = "input.conditionedPanels==4",
                     uiOutput("downLink4")
    ),
    conditionalPanel(condition = "input.conditionedPanels==5 || input.conditionedPanels==4",
                     sliderInput("ko_sd", "SD cutoff for KO terms", value = 20, min = 0, max = 100,step = 0.1)
    ),
    conditionalPanel(condition = "input.conditionedPanels==5",
                     actionButton("down5", "Download KEGG map")
    ),
    conditionalPanel(condition = "input.conditionedPanels==6",
                     fileInput('InFile', 'Upload previously saved Rdata file.'),
                     actionButton("loadRdata", "Upload Rdata"),
                     textInput("Rdataname","Enter file name for Rdata being saved (with '.Rdata' in the end"),
                     actionButton("saveRdata", "Save current Rdata")
    ),
    conditionalPanel(condition = "input.conditionedPanels==7",
                     h3("Your Metadata"),
                     #selectInput("colName", ColNameSelectorTitle, choices = colnames(mdt), selected = colName, selectize = FALSE),
                     uiOutput("colNameO")
    ),
    width = 3),

  mainPanel(
    tabsetPanel(
      tabPanel(titletab1, uiOutput("dynamic1") ,uiOutput("dynamic11"),uiOutput("dynamic111"), value = 1),
      tabPanel(titletab2, uiOutput("dynamic2"),uiOutput("dynamic21"),uiOutput("dynamic211"), value = 2),
      tabPanel(titletab3, uiOutput("dynamic3"),uiOutput("dynamic31"),uiOutput("dynamic311"), value = 3),
      tabPanel(titletab4, uiOutput("dynamic4"),uiOutput("dynamic41"),uiOutput("dynamic411"), value = 4),
      tabPanel(titletab5, imageOutput("Pathway",width = "100%", height = "400px"), value = 5),
      id = "conditionedPanels",
      tabPanel(titletab6, selectInput(inputId = "set_taxlevel1" , label = taxone, choices = tax1n, selected = tax1selected),
               selectInput(inputId = "set_taxlevel2", label = taxtwo, choices = tax2n, selected = tax2selected),
               selectInput(inputId = "set_funlevel1", label = funcone, choices = func1n, selected = func1selected),
               selectInput(inputId = "set_funlevel2", label = functwo, choices = func2n, selected = func2selected),
               selectInput(inputId = "colorPalette", label = labelColPal, choices = rownames(brewer.pal.info), selected = currentPalette),
               actionButton("showAllCols", "Show All Color Palettes"),
               plotOutput("paletteOutput"),
               actionButton("save_changes", "Save Changes"), value = 6),
      tabPanel(titletab7,
               DT::DTOutput("metadata"),
               # rHandsontableOutput("hot"),
               # div(class='row', div(h3("Edit Your Metadata"),
               #                      class="col-sm-5", uiOutput("ui_newcolname")),
               #     div(class="col-sm-4", h3("Select the type of a new column"),
               #         radioButtons("newcolumntype", "Type of a new column",
               #                      c("integer", "double", "character"))),
               #     div(class="col-sm-3")
               # ),
               # actionButton("addcolumn", "Create a new column"),
               # actionButton("addcolumn2", "Save a new column"),
               # actionButton("saveBtn", "Save"),
               value = 7) #dataTableOutput("table1")
    ), width = 9)
))
