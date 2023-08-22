editableUI <- function(id,
                       width=800,height=600) {
  ns <- NS(id)
  cat('editableUI',id,'\n')
  cat('editableUI',class(ns("res")),'\n')
  tagList(
    textInput(ns("colName"), "Column name"),
    actionButton(ns("create"), "Create Column", class = "pull-left btn-info"),
    #actionButton(ns("save"), "Save table", class = "pull-right btn-info"),
    downloadButton(ns('saveBtn'),'Save dataset'),
    #shinySaveButton(ns("save"), "Save file", "Save file as ...", filetype=list(rdata="Rdata"), class = "pull-right btn-info"),
    DT::DTOutput(ns("res"),height = '100%')
    #uiOutput("dynamic4")
    # rHandsontableOutput((rhandsontable(ns("res"), width = 600, height = 300,readOnly = TRUE) %>%
    #                        hot_cols(fixedColumnsLeft = 1) %>%
    #                        hot_rows(fixedRowsTop = 1)),width = 600)
  )
}

getTableDef<-function(pool,mgList){
  con<-poolCheckout(pool)
  #as.data.frame(table,stringsAsFactors=FALSE)
  res<-as.data.frame(con %>% dbReadTable('smpl') %>%
                       dplyr::filter(mgrastid %in% mgList) %>%
                       dplyr::select(c('name','mgrastid','description')),
                     stringsAsFactors=FALSE)
  res[which(is.na(res),arr.ind = TRUE)]<-''
  poolReturn(con)
  return(res)
}
editable <- function(input, output, session, mgList,makeEmptyCol,
                     updateTable,getTable=getTableDef,
                     colWidths = c(350,150,500),
                     width=800,height=600) {

  # observeEvent(tbls(), {
  #   updateSelectInput(session, "tableName", choices = tbls())
  # })
  #
  values <- reactiveValues()

  output$saveBtn <- downloadHandler(
    filename = function() {
      "pathview.Rdata"
    },
    content = function(file) {
      volumes <- c("UserFolder"=getwd())
      con<-poolCheckout(pool)
      savePathview(con = con,mdt = values$resDT,filename = file)
      showModal(modalDialog(
        title = "Your data have been saved",
        "Congratulations, Your data has been saved",
        easyClose = TRUE, footer = NULL
      ))
      poolReturn(con)
    }
  )

  # observeEvent(input$save, {
  #   volumes <- c("UserFolder"=getwd())
  #   shinyFileSave(input, "save", roots=volumes, session=session)
  #   fileinfo <- parseSavePath(volumes, input$save)
  #   con<-poolCheckout(pool)
  #   savePathview(con = con,mdt = values$resDT,filename = fileinfo$datapath)
  #   showModal(modalDialog(
  #     title = "Your data have been saved",
  #     "Congratulations, Your data has been saved",
  #     easyClose = TRUE, footer = NULL
  #   ))
  #   poolReturn(con)
  # })
  # # observe({
  # #   req(tabName)
  # #   #    req(tabName %in% db_list_tables(pool))
  # #   # cols <- db_query_fields(pool, input$tableName)
  # #   # updateCheckboxGroupInput(session, "select",
  # #   #                          choices = cols, selected = cols, inline = TRUE)
  # # })

  observeEvent(input$create, {
    if (nchar(input$colName) < 1) {
      showModal(modalDialog(
        title = "No column specified",
        "Please provide the name of column",
        easyClose = TRUE, footer = NULL
      ))
      return()
    }
    if(input$colName %in% names( values$resDT)){
      showModal(modalDialog(
        title = "Such column exists",
        "Please provide new name for the column",
        easyClose = TRUE, footer = NULL
      ))
      return()
    }
    cat("observeEvent(input$create",dim(values$resDT),'\n')
    if(is.null(values$saved)) values$saved=FALSE
    if(!values$saved){
      values$oldData<-values$resDT
    }
    values$resDT[,input$colName]<-NA
    cat("observeEvent(input$create",dim(values$resDT),'\n')
    values$saved=FALSE
  })

  # observeEvent(input$save, {
  #   #req(checkTable(values$resDF))
  #  # con<-poolCheckout(pool)
  #   cat('input$save\n')
  #   need(tryCatch(updateTable(values$oldData,values$newData)),"table update was unsuccessful")
  #   #values$resDT<-getTable(con,mgList)
  #   values$oldData<-values$resDT
  #   #poolReturn(con)
  #   values$saved=TRUE
  # })
  #
  dt<-reactive({
    values$resDT
  })
  # observe({
  #   reqTable(input$tableName)
  #   req(input$select)
  #   updateSelectInput(session, "filter", choices = input$select)
  # })
  observe({
    if(!is.null(input$res)){
      cat('class(hot_to_r',class(values$resDT),'\n')
      values$newData <- values$resDT#hot_to_r(input$res)
      values$saved=FALSE
    }
  })

  # observe({
  #   reqTable(input$tableName)
  #   req(input$select)
  #   updateSelectInput(session, "filter", choices = input$select)
  # })
  #
  # observe({
  #   reqColInTable(input$tableName, input$filter)
  #   df <- as_data_frame(pool %>% tbl(tabName))
  #   allUniqueVals <- unique(df[[input$filter]])
  #   updateCheckboxGroupInput(session, "vals",
  #                            choices = allUniqueVals, selected = allUniqueVals, inline = TRUE)
  # })

  proxy = DT::dataTableProxy('res')

  observeEvent(input$res_cell_edit, {
    values$saved=FALSE
    cat('input$res_cell_edit\n')
    info = input$res_cell_edit
    str(info)
    i = info$row
    j = info$col
    v = info$value
    values$resDT[i, j] <<- DT::coerceValue(v, values$resDT[i, j])
    values$newData<-values$resDT
    cat('toProxy',pander(head(values$newData),caption = "Proxy"),'\n')
    DT::replaceData(proxy, values$resDT, resetPaging = FALSE)  # important
  })

  #output$res <- renderRHandsontable({
  output$res <- DT::renderDT({
    #  values$resDT
    #   reqColInTable(input$tableName, input$filter)
    if(is.null(values$resDT)){
      #con<-poolCheckout(pool)
      cat('output$res: tblExsists',mgList,'\n')

      # filterVar <- sym(input$filter)
      # vals <- input$vals
      values$resDT<-getTable(pool,mgList)
      values$oldData<-values$resDT
      values$saved=TRUE
      #poolReturn(con)
    }else if(dim(values$resDT)[1]<length(mgList)){
      cat('add new row\n')
      idx<-match(mgList,values$resDT$mgrastid)
      add<-getTable(con,mgList[is.na(idx)])
      if(dim(add)[2]!=dim(values$resDT$mgrastid)[2]){
        idx<-match(names(values$resDT),names(add))
        add[names(values$resDT)[is.na(idx)]]<-NA
      }
      values$resDT<<-rbind(values$resDT,add)
    }
    values$resDT<<-unique(values$resDT)
    cat('editable output$res',dim(values$resDT),length(mgList),'\n')
    cat('editable output$res',class(values$resDT),'\n')
    #mat<-as.matrix(values$resDT)
    #cat('toMatrix',pander(head(values$resDT),caption = "Matrix"),'\n')
    #values$rhRes<-rhandsontable(as.matrix(values$resDT), width = width, height = height,search=TRUE,readOnly = FALSE,useTypes = FALSE) %>%
#     values$rhRes<-rhandsontable(values$resDT, width = width, height = height,search=TRUE,useTypes = FALSE) %>%
#       hot_col("mgrastid", readOnly = TRUE)%>%
#       hot_col(1:dim(values$resDT)[2], strict=FALSE,allowInvalid=TRUE) %>%
#       hot_cols(colWidths = colWidths) %>%
#       hot_cols(fixedColumnsLeft = 1) %>%
# #      hot_rows(fixedRowsTop = 1)%>%
#       #hot_context_menu(allowRowEdit = TRUE, allowColEdit = TRUE) #%>%
#       hot_cols(columnSorting = TRUE)
#     return(values$rhRes)
    return(datatable(values$resDT,editable = TRUE))
  })
}


savePathview<-function(con,mdt,filename='pathview.RData'){
 rownames(mdt)<-mdt$mgrastid
 funtaxall<-prepFunTaxAll(con,mdt)
 idx<-match(getSampleNames(funtaxall),rownames(mdt))
 cat('loadData-1:',idx,'\n')
 mdt<-mdt[idx,]
 d.kres<-prepKEGG(con,mdt)
 cat('loadData-2:\t',rownames(mdt),'\n\t\t',getSampleNames(funtaxall),'\n')

 cat('save',getwd(),filename,'\n')
 load("../../data/keggmappings.rda")
 save(funtaxall, d.kres, mdt, kegg, ko.path.name, file = filename)
}

prepFunTaxAll<-function(conn,mdt){
  cat(format(Sys.time(), "%b %d %X"),'\tStart SEED\n')
  res<-list()
  system.time(
    for(m in mdt$mgrastid){
      cat(format(Sys.time(), "%b %d %X"),'\tStart:',m,'\n')
      smpl<-dbGetQuery(conn,'select distinct md5sum,mgrastid,abundance,usp,ufun from mv_funtaxab where mgrastid=?',m);
      d.res<-data.table(smpl)
      names(d.res)[3]<-'ab'
      cat(format(Sys.time(), "%b %d %X"),'\tprepFTA',m,dim(d.res),length(res),dim(mdt),'\n')
      d.ab<-getAbundanceMD5FromDT(d.res)
      d.ab$mgid<-m
      res[[length(res)+1]]<-d.ab
      cat(format(Sys.time(), "%b %d %X"),'\tDone:',m,'\n')
    }
  )
  resDF<-do.call(rbind,res)
  d.bm <- our.aggregate(resDF)
  usp<-unique(d.bm$usp)
  taxL<-dbGetQuery(conn,'select taxid,usp,name,species,genus,family,ordr,class,phylum,dom from taxlite')
  names(taxL)<-c("taxid","usp","name","species","genus","family","order","class","phylum","domain")
  taxL<-taxL[match(usp,taxL$usp),]
  taxall<-merge(d.bm,data.table(taxL),by.y = 'usp',by.x='usp')
  funtree<-dbGetQuery(conn,'select accession,fun1,fun2,fun3,fun4 from seedlite')
  names(funtree)<-c('id',"FUN4",'FUN3','FUN2','FUN1')
  funtaxallL <- merge(taxall, funtree, by.x = 'ufun', by.y = 'id')
  funtaxall<-funtaxallL[,c("usp", "species","genus","family","order","class","phylum","domain","md5","FUN1","FUN2","FUN3","FUN4",grep('mgm',names(funtaxallL),value = TRUE)), with=FALSE]
  names(funtaxall)[10]<-'ufun'
  return(funtaxall)
}

prepKEGG<-function(con,mdt){
  res<-list()
  cat(format(Sys.time(), "%b %d %X"),'\tStart KEGG\n')
  system.time(
    for(m in mdt$mgrastid){
      cat(format(Sys.time(), "%b %d %X"),'\tStart:',m,'\n')
      smpl<-dbGetQuery(con,'select distinct id,md5,annotation,ko from mv_dkres where id=?',m);
      d.res<-data.table(smpl)
      res[[length(res)+1]]<-d.res
      cat(format(Sys.time(), "%b %d %X"),'\tDone:',m,'\n')
    }
  )
  resDF<-do.call(rbind,res)
  names(resDF)<-c('.id', "md5","annotation","ko")
  return(resDF)
}