#
# Shiny Nanostring25
# images components
#


output$boxplotIcon <- renderImage({
    filename <- normalizePath(file.path('pics','boxplotIcon.png'))          
    # Return a list containing the filename and alt text
    return(list(src = filename,
                alt = "boxplot view"))       
}, deleteFile = FALSE)

output$refValuesIcon <- renderImage({
    filename <- normalizePath(file.path('pics','refValuesIcon.png'))          
    # Return a list containing the filename and alt text
    return(list(src = filename,
                alt = "Referance values table"))       
}, deleteFile = FALSE)

output$pcaIcon <- renderImage({
    filename <- normalizePath(file.path('pics','pcaIcon2.png'))          
    # Return a list containing the filename and alt text
    return(list(src = filename,
                alt = "PCA view"))       
}, deleteFile = FALSE)

output$vennIcon <- renderImage({
    filename <- normalizePath(file.path('pics','vennIcon.png'))          
    # Return a list containing the filename and alt text
    return(list(src = filename,
                alt = "Venn view"))       
}, deleteFile = FALSE)