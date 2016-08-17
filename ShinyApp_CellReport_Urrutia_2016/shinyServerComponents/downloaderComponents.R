#
# Shiny Nanostring25
# Downloader Components
#

#
# mRNA Boxplot Download
#
output$downloadmRNABoxplot <- downloadHandler(
    
    # Filename 
    filename = function() {
        paste0("LabExMI_Nanostring25_mRNABoxplots_", paste0(input$boxplot_mRNAs, collapse="-"), "_", format(Sys.Date(),"%Y-%m-%d"), ".png")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
        p <- facettedmRNABoxplots(df.data = df.data.long.filtered(), df.stimulus.color = df.stimulus.color, compactPlotFlag = input$boxplot_muliple_labels)
        # Write to a file specified by the 'file' argument
        ggsave(plot=p, filename=file,width = 15, height=getBoxPlotHeight()/80)
    }
)

# downloadHandler() takes two arguments, both functions.
# The content function is passed a filename as an argument, and
#   it should write out data to that filename.
output$downloadReferenceValues <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
        paste0("LabExMI_Nanostring25_ReferenceValues_",input$valueType,"_", format(Sys.Date(),"%Y-%m-%d"), ".tsv")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
        #sep <- switch(input$filetype, "csv" = ",", "tsv" = "\t")
        
        # Write to a file specified by the 'file' argument
        write.table(df.ref.values.selected(), file, sep = "\t",
                    row.names = FALSE)
    }
)