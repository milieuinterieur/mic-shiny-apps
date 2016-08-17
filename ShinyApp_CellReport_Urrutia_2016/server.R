# 
# Author: Vincent Rouilly
# Oct 2015
#
# SERVER.R 
#

#
# shiny server code
#
shinyServer(function(input, output){
    
    ##############################################
    # load reactive expressions
    source("shinyServerComponents/reactives.R", local = TRUE)
    
    ##############################################
    # create UI input widgets from data
    source("shinyServerComponents/uiWidgetComponents.R", local = TRUE)
    
    ##############################################
    # load plot components
    source("shinyServerComponents/plotComponents.R", local = TRUE)
    
    ##############################################
    # load data table components
    source("shinyServerComponents/dataTableComponents.R", local = TRUE)
   
    ##############################################
    # load downloader components
    source("shinyServerComponents/downloaderComponents.R", local = TRUE)
    
    ##############################################
    # load images components
    source("shinyServerComponents/imageComponents.R", local = TRUE)
    
})
