dtTable <- function(dataFromDB, editable) {
        DT::datatable(
                data = dataFromDB,
                rownames = FALSE,
                
                extensions = list(
                        # Bottons
                        'Buttons'= NULL,
                        
                        # The Scroller extension makes it possible to only render
                        # the visible portion of the table.
                        'Scroller' = NULL
                ),
                
                options = list(
                        columnDefs = list(list(
                                targets = c(0:(length(colnames(dataFromDB)) - 1)),
                                render = JS(
                                        "function(data, type, row, meta) {",
                                                "return type === 'display' && !!data && data.length >= 15 ?",
                                                "'<span title=\"' + data + '\">' + data.substr(0, 14) + '..</span>' : data;",
                                        "}"
                                )
                        )),
                        
                        # Bottons
                        dom = 'Bfrtip',
                        buttons = list(
                                'colvis',
                                list(
                                        extend = 'collection',
                                        buttons = c('csv', 'excel', 'pdf'),
                                        text = 'Download'
                                )
                        ),
                        
                        # The Scroller extension makes it possible to only render
                        # the visible portion of the table.
                        deferRender = TRUE,
                        scrollX = TRUE,
                        scrollCollapse = TRUE,
                        paging = FALSE,
                        scrollY = "61vh"
                ),
                
                selection = "none",
                editable = editable
        )
}
