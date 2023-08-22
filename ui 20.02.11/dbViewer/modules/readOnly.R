readOnlyUI <- function(id,
                       width=800,height=800) {
  ns <- NS(id)
  cat('readOnlyUI',id,'\n')
  cat('readOnlyUI',class(ns("res")),'\n')
  tagList(
    # tableOutput(ns("res"))
    rHandsontableOutput(ns("res"),width = width,height=height), br(),
    verbatimTextOutput(ns('resultDisplay1')),
    verbatimTextOutput(ns('resultDisplay2')),
    verbatimTextOutput(ns('resultDisplay3'))

    # rHandsontableOutput((rhandsontable(ns("res"), width = 600, height = 300,readOnly = TRUE) %>%
    #                        hot_cols(fixedColumnsLeft = 1) %>%
    #                        hot_rows(fixedRowsTop = 1)),width = 600)
  )
}

readOnly <- function(input, output, session, pool, tabName,
                     colWidths = c(50,150,300),width=800,height=600) {

  # observeEvent(tbls(), {
  #   updateSelectInput(session, "tableName", choices = tbls())
  # })
  #

  observe({
    req(tabName)
    req(tabName %in% dbListTables(pool))
    # cols <- db_query_fields(pool, input$tableName)
    # updateCheckboxGroupInput(session, "select",
    #                          choices = cols, selected = cols, inline = TRUE)
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

  output$res <- renderRHandsontable({
 #   reqColInTable(input$tableName, input$filter)
    cat('output$res: tblExsists',tabName,(tabName %in% dbListTables(pool)),'\n')

    # filterVar <- sym(input$filter)
    # vals <- input$vals

    res<-dbReadTable(pool,tabName)
    cat('output$res',dim(res),'\n')
    cat('output$res',class(res),'\n')
    rhRes<-rhandsontable(res, width = width, height = height,search=TRUE,selectCallback = TRUE,readOnly = TRUE) %>%
      hot_cols(colWidths = colWidths) %>%
      hot_cols(fixedColumnsLeft = 1) %>%
      hot_rows(fixedRowsTop = 0)%>%
      hot_cols(columnSorting = TRUE)
    return(rhRes)

  })
  # output$resultDisplay1 <- renderPrint({
  #   cat('Selected Row      : ', input$res_select$select$r,'\n')
  #   cat('Selected Column   : ', input$res_select$select$c,'\n')
  #   cat('Selected Cell Val : ', input$res_select$data[[input$res_select$select$r]][[input$res_select$select$c]],'\n')
  #   cat('Selected Range    : ','(',input$res_select$select$r, ',',input$res_select$select$c, ')~(',
  #       input$res_select$select$r2,',',input$res_select$select$c2,')', sep="",'\n')
  # })
  # # Display Changed Row
  # output$resultDisplay2 <- renderPrint({
  #   cat('Changed Cell Row Column : ', input$res$changes$changes[[1]][[1]], input$res$changes$changes[[1]][[2]],'\n')
  #   cat('Changed Cell Old Value  : ', input$res$changes$changes[[1]][[3]],'\n')
  #   cat('Changed Cell New Value  : ', input$res$changes$changes[[1]][[4]])
  # })
  #
  # res_select_r <- reactive(input$res_select$select)
  # observe({
  #     selectedRow <- res_select_r()
  #     idx <- selectedRow$r
  #     if(is.null(idx)) {
  #         env4$select_idx <<- NULL ; cat("====> NULL")
  #     } else {
  #         env4$select_idx <<- idx  ; cat("====>",env4$select_idx)
  #     }
  # })

  # output$resultDisplay3 <- renderPrint({
  #   env4$select_idx <<- input$res_select$select$r
  #   cat('Selected Row : ', env4$select_idx,'\n')
  # })

  select<-reactive({input$res_select$data[[input$res_select$select$r]][[5]]})
#  select<-reactive({input$res_select$data[[input$res_select$select$r]][[4]]})
  return(select)
}

# init_env4 <- function(){
#   env <- new.env()
#   env$select_idx <- NULL
#   return(env)
# }
# env4 <- init_env4()
#
