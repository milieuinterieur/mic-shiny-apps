#
# Shiny Nanostring25
# Data Table Components
#


#
# Render Reference Values Table
#
output$dataTable <- DT::renderDataTable({
    req(df.ref.values.selected, input$ref_table_stimuli, input$refTableMaxQvalue)
    
    df.data.selected <- df.ref.values.selected()
    
    # stick q value data for highlights
    startIndex <- dim(df.data.selected)[2] + 1
    endIndex <- startIndex + dim(df.data.selected)[2] - 2
    
    df.data.selected.highlight <- df.ref.table.qvalue[,c(input$ref_table_stimuli)]
    names(df.data.selected.highlight) <- paste0(names(df.data.selected.highlight),".highlight")
    df.data.selected <- cbind(df.data.selected, df.data.selected.highlight)
    DT::datatable(df.data.selected, filter = "top", options = list(
        columnDefs = list(list(targets = startIndex:endIndex , visible = FALSE)),
        pageLength = 25)) %>%
        formatStyle(
            input$ref_table_stimuli, names(df.data.selected.highlight),
            color = styleInterval(input$refTableMaxQvalue, c('red', 'black')),
            backgroundColor = styleInterval(input$refTableMaxQvalue, c('lightblue', 'white'))
        )
})

#
# Venn Diagram Data Table
#
output$vennDiagramDT <- DT::renderDataTable({
    getVennDiagramDataTable(venn.id = input$venn_id)  
}, filter = "top")