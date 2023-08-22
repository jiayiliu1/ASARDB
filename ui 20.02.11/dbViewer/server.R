#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(RColorBrewer)
library(ggplot2)
#library(DBI)
#library(dplyr)
#library(dtplyr)
library(tibble)
library(pool)
library(rlang)
library(MonetDB.R)
library(rhandsontable)
library(pander)
library(data.table)
library(asarDB)
library(shiny)
library(shinyFiles)

#source('db.R')
source("./modules/readOnly.R", local = TRUE)
source("./modules/editable.R", local = TRUE)

options(shiny.maxRequestSize=500*1024^2)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  mgList<-c()
#  callModule(readOnly, "prj", pool,tabName="study",colWidths=c(50,100,100,500,50))
  selM<-callModule(readOnly, "mgen", pool,tabName="mv_smpl_stat",colWidths=c(250,100,50,300,100,100,100),width=1200)
  #selM<-callModule(readOnly, "mgen", pool,tabName="smpl",colWidths=c(50,200,50,100,300,100,100))
  observe({
    cat('SelM',class(selM()),'"',selM(),'"\n',mgList,'\n')
    mgList<<-unique(c(mgList,selM()))
    if(length(mgList)>0){
    output$selected<-renderPrint(paste(mgList,collapse = ','))
    }
    cat('SelM',class(selM()),'"',selM(),'"\n')
  })
  callModule(editable, "dtset",mgList=mgList, colWidths=c(50,200,500),
             makeEmptyCol=makeDset,updateTable=updateDset)
  callModule(readOnly, "stat", pool,tabName="smpl",colWidths=c(50,200,50,100,300,100,100))
  #callModule(readOnly, "stat", pool,tabName="mv_smpl_stat",colWidths=c(100,100,50,300,100,100,100))
})

updateDset<-function(olddata,newdata){
  cat('Update dataset\n',pander(head(olddata),caption = "old"),'\n=====\n',pander(head(newdata),caption = "new"),'\n====\n')
}

makeDset<-function(){
  cat('makeDset\n')
}