#
# Shiny Nanostring25
# reactives.R
#
#

#
# Figure PCA Settings
#

#default PCA projection settings
pca.settings <- reactiveValues(space = "I3T", projected = c('IFN-A'), genes = "572genes")

# Fig1A settings
observeEvent(input$pca.fig1A, {
    pca.settings$space <- "I3T"
    pca.settings$projected <- c()
    pca.settings$genes <- "572genes"
})

# Fig1B settings
observeEvent(input$pca.fig1B, {
    pca.settings$space <- "I3T"
    pca.settings$projected <- c()
    pca.settings$genes <- "44genes"
})

# FigS3A settings
observeEvent(input$pca.figS3A, {
    pca.settings$space <- "I3T"
    pca.settings$projected <- c('IFN-A')
    pca.settings$genes <- "44genes"
})

# FigS5A settings
observeEvent(input$pca.figS4A, {
    pca.settings$space <- "I3T"
    pca.settings$projected <- c('FSL','pIC','LPS','FLA','GARD','r848','ODN')
    pca.settings$genes <- "44genes"
})

# FigS5A settings
observeEvent(input$pca.figS4B, {
    pca.settings$space <- "I3T"
    pca.settings$projected <- c('WGP','LAM','CPPD')
    pca.settings$genes <- "44genes"
})

# # Fig4A settings
# observeEvent(input$pca.fig4A, {
#     pca.settings$space <- "I3T"
#     pca.settings$projected <- c('FSL','pIC','LPS','FLA','GARD','r848','ODN')
#     pca.settings$genes <- "44genes"
# })
# 
# # Fig4B settings
# observeEvent(input$pca.fig4B, {
#     pca.settings$space <- "I3T"
#     pca.settings$projected <- c('FSL','pIC','LPS','FLA','GARD','r848','ODN')
#     pca.settings$genes <- "44genes"
# })

# FigS7 settings
observeEvent(input$pca.figS5, {
    pca.settings$space <- "I3T"
    pca.settings$projected <- BUG.stimuli
    pca.settings$genes <- "44genes"
})


#
# Figure Boxplot Settings
#

boxplot.settings <- reactiveValues(stimuli = c('Null',I3T.stimuli), mRNAs = c('CCL8','CXCL9','IL6','CD83'))

# FigSXXX settings
observeEvent(input$boxplot.figSXXX, {
    boxplot.settings$stimuli <- c('Null',I3T.stimuli)
    boxplot.settings$mRNAs <- top44.genes[1:5]
})

# Fig2A settings
observeEvent(input$boxplot.fig2A, {
    boxplot.settings$stimuli <- c('Null',I3T.stimuli)
    boxplot.settings$mRNAs <- c('IFNB1','IFNG','IL1B','TNF')
})

# Fig2C settings
observeEvent(input$boxplot.fig2C, {
    boxplot.settings$stimuli <- c('Null',I3T.stimuli)
    boxplot.settings$mRNAs <- c('CCL8','CXCL9','IL6','CD83')
})

# Fig3B settings
observeEvent(input$boxplot.fig3B, {
    boxplot.settings$stimuli <- c('Null',TLR.stimuli)
    boxplot.settings$mRNAs <- c('IFNB1','IFNG','IL1B','TNF')
})

#
# Reactive Reference Values Datasets
#
df.ref.values.selected <- reactive({ 
    if(is.null(input$valueType))
        return
    
    if(input$valueType =="qvalue"){
        df.data <- getQValues()
        df.data.selected <- df.data[,c('variable',input$ref_table_stimuli)]
    } else if(input$valueType =="fc"){
        df.data <- getFCValues()
        df.data.selected <- df.data[,c('variable',input$ref_table_stimuli)]
    } else if(input$valueType =="median"){
        df.data <- getMedianValues()
        df.data.selected <- df.data[,c('variable',input$ref_table_stimuli)]
    } else if(input$valueType =="cv"){
        df.data <- getCvValues()
        df.data.selected <- df.data[,c('variable',input$ref_table_stimuli)]
    }
    return(df.data.selected)
}) 

#
# Filter dataset based on Stimuli and mRNAS
#
df.data.long.filtered <- reactive({
    df.data.long.filtered <- df.data.long %>% 
        filter(StimulusName %in% input$boxplot_stimuli, 
               mRNA %in% input$boxplot_mRNAs)
    df.data.long.filtered$mRNA <- factor(as.character(df.data.long.filtered$mRNA), levels=input$boxplot_mRNAs)
    df.data.long.filtered$StimulusName <- factor(as.character(df.data.long.filtered$StimulusName), levels=input$boxplot_stimuli)
    
    df.data.long.filtered                             
})

#
# Compute PCA components
#
pca.meta.info <- c('StimulusName','Color','Donor.ID')

pca.fond.stimuli <- reactive({
    req(input$pca_stimuli_space)
    
    pca_fond_stimuli <- switch(input$pca_stimuli_space,
                               I3T = I3T.stimuli,
                               TLRs = TLR.stimuli,
                               Microbes = BUG.stimuli)
#     if (input$pca_stimuli_space=='I3T'){
#         pca_fond_stimuli = I3T.stimuli
#     } else {
#         pca_fond_stimuli = TLR.stimuli
#     }
    pca_fond_stimuli
})

#
# PCA results
#
df.pca.coords <- reactive({    
    req(input$pca_gene_list, pca.fond.stimuli)
    getCoordinatesFromPcaAndProjections(nbDim=3,
                                        df.data = df.data.wide.DonorRegressed[,c(pca.meta.info, input$pca_gene_list)],
                                        pca.stimuli = pca.fond.stimuli(), 
                                        projected.stimuli = input$pca_proj_stimuli)
})

df.data.mds.ready.filtered.normalized <- reactive({
    
    #print(input$mrna_signature_stimuli)
    if(!is.null(input$mrna_signature_stimuli)) {
        df.data.mds.ready.filtered <- subset(df.data.mds.ready, select= input$mrna_signature_stimuli)
    } else {
        df.data.mds.ready.filtered <- subset(df.data.mds.ready, select= c('LPS','FLA'))
    }
    #print(head(df.data.mds.ready.filtered))
    # helper function
    if (input$metricNorm=="MaxNorm"){
        max.abs <- function(x) max(abs(x))
        # extract absolute max value
        rowmax.values <- apply(df.data.mds.ready.filtered, 1, max.abs)
        # normalize by absolute max value
        df.data.mds.ready.filtered.normalized2 <- sweep(df.data.mds.ready.filtered, 1, rowmax.values, "/")
        #print(head(df.data.mds.ready.filtered.normalized2))
    } else {
        sd.values <- apply(df.data.mds.ready.filtered, 1, sd)
        df.data.mds.ready.filtered.normalized2 <- sweep(df.data.mds.ready.filtered, 1, sd.values, "/")
    }
    
    return(df.data.mds.ready.filtered.normalized2)
})

mds.model <- reactive({
    # compute distance
    d <- dist(df.data.mds.ready.filtered.normalized())
    # compute MDS fit
    fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
    return(fit)
})

df.mds.points <- reactive({
    df.data.mds <-  data.frame(coordX = mds.model()$points[,1],
                               coordY = mds.model()$points[,2],
                               names = row.names(df.data.mds.ready.filtered.normalized()))
    return(df.data.mds)
})

selected.mRNA.signatures <- reactive({
    brushedPoints(df.mds.points(), input$plot_brush, 
                   xvar ="coordX", 
                   yvar= "coordY",
                   panelvar1 = NULL, 
                   panelvar2 = NULL, 
                   allRows = FALSE)
})

nb.boxplot.rows <- reactive({
    nb.rows <- (length(selected.mRNA.signatures$names)/5)*(200)
    return(nb.rows)
})
