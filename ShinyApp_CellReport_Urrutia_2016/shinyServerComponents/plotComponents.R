#
# Shiny Nanostring25
# Plot Components
#


# helper function to deal with plot height setting
getBoxPlotHeight <- function(){
    return(length(input$boxplot_mRNAs)*(180+ifelse(input$boxplot_muliple_labels, 0,60)))
}


#
# mRNAs boxplot view with selected Stimuli
#
output$mRNAsBoxplots <- renderPlot({        
    
    p <- facettedmRNABoxplots(df.data = df.data.long.filtered(), df.stimulus.color = df.stimulus.color, compactPlotFlag = input$boxplot_muliple_labels)
    print(p)
}, height=getBoxPlotHeight)

#
# Venn Diagram Plot
#
output$vennDiagramPlot <- renderPlot({
    req(input$venn_id, input$vennMaxQvalue)
    #plotVennDiagram(input$venn_id)
    plotDynamicVennDiagram(input$venn_id, maxQValue = input$vennMaxQvalue)
})


# helper function to deal with plot height setting
getPcaPlotHeight <- function(){
    return(1*(350))
}
#
# PCA View
#
output$pcaPlots <- renderPlot({
    req(df.pca.coords, pca.fond.stimuli)
    
    df.pca.3dcoords <- df.pca.coords()
    # re-order Stimulus Names
    df.pca.3dcoords$StimulusName <- factor(as.character(df.pca.3dcoords$StimulusName), levels = c(pca.fond.stimuli(), input$pca_proj_stimuli))
    
    #color palette
    color.palette <- brewer.pal(n = 12, name='Paired')
    
    #print(paste0("pca_proj_stimuli: ", input$pca_proj_stimuli))
    #print(paste0("pca_fond_stimuli: ", pca.fond.stimuli()))
    # Pile on PC1vsPC2 - 192 genes - with labels
    pca.plot <- facettedProjectedPcaPlot(df.data=df.pca.3dcoords[,c('StimulusName','Donor.ID','PC1','PC2')], 
                                         xName='PC1', yName='PC2',
                                         ellipse.stimuli=c(pca.fond.stimuli(),input$pca_proj_stimuli), #input.stimuli, 
                                         points.stimuli=c(),
                                         stimuli.colors=c(color.palette[1:length(pca.fond.stimuli())], rep('magenta',times = ifelse(is.null(input$pca_proj_stimuli),0,length(input$pca_proj_stimuli)))),
                                         #stimuli.colors=c('olivedrab3','purple','lightsteelblue3','turquoise1','red','red','red','red', 'red', 'red', 'red'),
                                         labels=TRUE, ellipse.linetype=1, #c(rep(1,times=length(input$pca_fond_stimuli)), 
                                         #rep(2, times=length(input$pca_proj_stimuli))), 
                                         xy.limits=c(-15,15))
    print(pca.plot)
}, height=getPcaPlotHeight)

output$pcaPlots2 <- renderPlot({
    req(df.pca.coords, pca.fond.stimuli)
    
    df.pca.3dcoords <- df.pca.coords()
    # re-order Stimulus Names
    df.pca.3dcoords$StimulusName <- factor(as.character(df.pca.3dcoords$StimulusName), levels = c(pca.fond.stimuli(), input$pca_proj_stimuli))
    #df.scatterplot.selectedgenes$StimulusName <- factor(as.character(df.scatterplot.selectedgenes$StimulusName), levels = c(projected.stimuli, input.stimuli))
    
    #color palette
    color.palette <- brewer.pal(n = 12, name='Paired')
    
    # Pile on PC1vsPC2 - 192 genes - with labels
    pca.plot <- facettedProjectedPcaPlot(df.data=df.pca.3dcoords[,c('StimulusName','Donor.ID','PC1','PC3')], 
                                         xName='PC1', yName='PC3',
                                         ellipse.stimuli=c(pca.fond.stimuli(),input$pca_proj_stimuli), #input.stimuli, 
                                         points.stimuli=c(),
                                         stimuli.colors=c(color.palette[1:length(pca.fond.stimuli())], rep('magenta',times = ifelse(is.null(input$pca_proj_stimuli),0,length(input$pca_proj_stimuli)))),
                                         #stimuli.colors=c('olivedrab3','purple','lightsteelblue3','turquoise1','red','red','red','red', 'red', 'red', 'red'),
                                         labels=TRUE, ellipse.linetype=1, #c(rep(1,times=length(input$pca_fond_stimuli)), rep(2, times=length(input$pca_proj_stimuli))), 
                                         xy.limits=c(-15,15))
    return(pca.plot)
}, height=getPcaPlotHeight)

#
# MDS Scatterplot
#
output$mdsScatterplot  <- renderPlot({
    
    p <- getMDSScatterplot(df.mds.data=df.mds.points(),
                           highlighted.mRNAs = input$mrna_signature_mRNAs)
    return(p)
})

# helper function to deal with plot height setting
getmRNASignatureBoxplotHeight <- function(){
    nb.boxplot.rows()
}
#
# mRNA signature across stimuli as boxplot
#
output$mRNASignatureBoxplots <- renderPlot({
    if(!input$mds_show_transformed_data){
        mRNAs <- as.character(selected.mRNA.signatures()$names)
        p <- getmRNASignatureAcrossStimuli(df.data.long=df.data.long, 
                                           selected.mRNAs=mRNAs,
                                           selected.stimuli=input$mrna_signature_stimuli)
    } else {
        df.transformed.data <- cbind(mRNA=row.names(df.data.mds.ready.filtered.normalized()),
                                     df.data.mds.ready.filtered.normalized())
        mRNAs <- as.character(selected.mRNA.signatures()$names)
        df.transformed.data.long <- df.transformed.data %>% #filter(mRNA %in% mRNAs) %>%
            gather("StimulusName",
                   "value",
                   one_of(names(df.data.mds.ready.filtered.normalized())))
        p <- getmRNANormalizedSignatureAcrossStimuli(df.data.long=df.transformed.data.long, 
                                                     selected.mRNAs=mRNAs,
                                                     selected.stimuli=input$mrna_signature_stimuli)
    }
    return(p)
}, height = "auto")

#
# Stimuli Correlations Matrix
#

output$stimuliCorrelationMatrix <- renderPlot({
    req(input$mrna_stimuli_correlations)
    
    require(corrplot)
    stimuli.filtered <- as.character(df.ttests.results %>% filter(mRNA==input$mrna_stimuli_correlations & FClog2>2 & qValue<0.001) %>%  
                                         select(one_of("StimulusName")) %>% unlist(use.names = FALSE))
    
    if(length(stimuli.filtered)>0){
        stimuli.filtered <- str_replace(stimuli.filtered, pattern = "-",replacement = ".")
        df.data.selected <- df.stimuli.corr[df.stimuli.corr$mRNA==input$mrna_stimuli_correlations, c('Stimulus','Null',stimuli.filtered)]
        row.names(df.data.selected) <- str_replace(as.character(df.data.selected$Stimulus), pattern = "-",replacement = ".")
        df.data.selected <- df.data.selected[row.names(df.data.selected) %in% c('Null',stimuli.filtered), c('Null',stimuli.filtered)]
        
        corrplot(as.matrix(df.data.selected), method="circle", shade.col=NA, tl.col="black", diag=FALSE,
                 col=colorRampPalette(c("blue","grey","grey","grey","red"))(200),
                 order="hclust", hclust.method="ward.D", tl.srt=45, mar = c(0,0,1,0))
        
        corrplot(as.matrix(df.data.selected),add=TRUE, type="lower", method="number",
                 order="hclust", hclust.method="ward.D", diag=FALSE, tl.pos="n", cl.pos="n",
                 col=colorRampPalette(c("blue","grey","grey","grey","red"))(200),
                 mar = c(0,0,1,0), title=paste0("Spearman correlations between stimulus conditions for ",input$mrna_stimuli_correlations," expression"))
        
    } else {
        p <- ggplot(data=data.frame())
        p <- p + geom_point()
        p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), axis.line = element_blank(), axis.text = element_blank())
        p <- p + xlim(0, 10) + ylim(0, 10)
        return(p)
    }
},width = 850, height=850)

output$stimuliCorrelationDendogram <- renderPlot({
    req(input$mrna_stimuli_correlations)
    
    stimuli.filtered <- as.character(df.ttests.results %>% filter(mRNA==input$mrna_stimuli_correlations & FClog2>2 & qValue<0.001) %>%  
                                         select(one_of("StimulusName")) %>% unlist(use.names = FALSE))
    if(length(stimuli.filtered)>0){
        stimuli.filtered <- str_replace(stimuli.filtered, pattern = "-",replacement = ".")
        
        df.data.selected <- df.stimuli.corr[df.stimuli.corr$mRNA==input$mrna_stimuli_correlations, -(1:2)]
        row.names(df.data.selected) <- str_replace(df.stimuli.corr[df.stimuli.corr$mRNA==input$mrna_stimuli_correlations, c('Stimulus')], pattern = "-",replacement = ".")
        df.data.selected <- 1 - abs(df.data.selected)
        fit <- hclust(as.dist(df.data.selected[row.names(df.data.selected) %in% c('Null',stimuli.filtered), c('Null',stimuli.filtered)]), method="ward.D")
        plot(fit, main=paste0("Hierarchical clustering based on Spearman correlation for ",input$mrna_stimuli_correlations," expression"), sub="", xlab="", mai=c(0.,0.,0.,0.), mar=c(0.,0.,0.,0.))
    } else {
        p <- ggplot(data=data.frame())
        p <- p + geom_point()
        p <- p + xlim(0, 10) + ylim(0, 10)
        p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       panel.background = element_blank(), axis.line = element_blank(),axis.text = element_blank())
        p <- p + annotate("text", label="No data to display here,", x=5, y=9, size=10)
        p <- p + annotate("text", label="as no stimulus condition has a ", x=5, y=7, size=10)
        p <- p + annotate("text", label=paste0("significant gene expression for ", input$mrna_stimuli_correlations), x=5, y=5, size=10)
        p <- p + annotate("text", label="Please select another gene.", x=5, y=2, size=6)
        return(p)
    }
    
}, width=800, height=400)

