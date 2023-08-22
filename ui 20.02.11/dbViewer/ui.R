#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(shiny)
library(shinydashboard)
library(rhandsontable)
source("./modules/readOnly.R", local = TRUE)
source("./modules/editable.R", local = TRUE)

# Define UI for application that draws a histogram
dashboardPage(
  dashboardHeader(title = 'Projects'),
  dashboardSidebar(width = 350,
                   sidebarMenu(
                     textInput("searchField", "Search"),
#                     menuItem("Projects",tabName = 'Projects'),
                     menuItem("Metagenomes",tabName = 'Metagenomes'),
                     menuItem("Datasets",tabName = 'Datasets'),
verbatimTextOutput('selected'),
                     menuItem("Statistics",tabName = 'Statistics')
                   )),
  dashboardBody(
    tabItems(
#      tabItem('Projects',readOnlyUI("prj")),
      tabItem('Metagenomes',readOnlyUI("mgen")),
      tabItem('Datasets',editableUI("dtset")),
      tabItem('Statistics',readOnlyUI("stat"))#,
    )
   ),
  title = "ASAR DB Viewer"
)
