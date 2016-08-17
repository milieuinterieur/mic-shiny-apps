#
# Shiny Nanostring25
# UI Widgets Components
#


# multi select for Stimuli, boxplot view
output$boxplotStimuliSelectControl <- renderUI({
    names <- levels(df.data.long$StimulusName)
    selectizeInput("boxplot_stimuli",
                   label = "Select Stimuli",
                   choices = names, 
                   multiple=TRUE,
                   #selected = c('Null','FSL','pIC','LPS','FLA','GARD','r848','ODN'),
                   selected = boxplot.settings$stimuli,
                   options = list(placeholder = 'select a stimulus'))
})

# multi mRNA select Control for boxplot view
output$boxplotmRNAsSelectControl <- renderUI({
    names <- unique(df.data.long$mRNA)
    selectizeInput("boxplot_mRNAs",
                   label="Select mRNAs (maximum of 5 at once.)", 
                   choices=names, 
                   multiple=TRUE,
                   #selected = c('IFNB1','IFNG','IL1B','TNF'),
                   selected = boxplot.settings$mRNAs,
                   options = list(placeholder = 'select an mRNA', maxItems = 5))
})

# multi select for Stimuli, boxplot view
output$pcaFondationalStimuliSelectControl <- renderUI({
    req( pca.settings$space)
    
    names <- c('[IFN-B, IFN-G, IL-1B, TNFa] based'='I3T',
               '[FSL, pIC, LPS, FLA, GARD, r848, ODN] based'='TLRs')#,
               #'[HKHP,HKSA,HKLR,HKEC,BCG,HKCA,IAV,SeV] Microbes' = 'Microbes')
    selectInput("pca_stimuli_space", 
                label = "Select PCA Stimuli Space",
                choices = names, 
                multiple=FALSE,
                #selected = 1
                selected = pca.settings$space)
})

# multi select for Stimuli, PCA view
output$pcaProjectedStimuliSelectControl <- renderUI({
    req(pca.settings, pca.fond.stimuli)
    
    names <- unique(df.data.long$StimulusName)
    stimuli.choices <- c("",setdiff(names,pca.fond.stimuli()))
    stimuli.selected <- switch(pca.settings$space, # input$pca_stimuli_space,
           I3T = TLR.stimuli,
           TLRs = BUG.stimuli[1],
           Microbes = TLR.stimuli[1])
    
    #if(input$pca_stimuli_space=='I3T'){
    #    stimuli.selected <- TLR.stimuli
    #} else if {
    #    stimuli.selected <- BUG.stimuli[1]
    #}
    selectInput("pca_proj_stimuli", 
                label = "Select Projected Stimuli (displayed in purple)",
                choices = stimuli.choices, 
                multiple=TRUE,
                #selected = stimuli.selected
                selected = pca.settings$projected)
})

# PCA gene set control
output$pcaGeneSetControl <- renderUI({
    req(pca.settings)
    
    choices <- switch(pca.settings$space, #input$pca_stimuli_space,
                               I3T = c("44 genes" = "44genes", "572 genes" = "572genes"),
                               TLRs = c("44 genes" = "44genes",
                                        "572 genes" = "572genes"),
                               Microbes = c("572 genes" = "572genes"))
    
    selectInput("pcaGeneSet", "Preset gene sets",
                choices,
                #selected=c('44genes')
                selected = pca.settings$genes)
    
})

# list of genes to use in PCA
output$pcaGeneSetList <- renderUI({
    req(input$pcaGeneSet)
    
    if (input$pcaGeneSet=='40random'){
        selected.genes <- sample(x = selected.mRNAs, size = 40, replace=FALSE)
    } else if (input$pcaGeneSet=='44genes'){
        selected.genes <- top44.genes
    #} else if (input$pcaGeneSet=='192genes'){
    #    selected.genes <- I3T.genes
    #} else if (input$pcaGeneSet=='202genes'){
    #    selected.genes <- I3T.202.genes
    } else {
        selected.genes <- selected.mRNAs
    }  
    selectInput("pca_gene_list",
                label = "Gene set used for PCA computation",
                choices = selected.mRNAs, 
                multiple=TRUE,
                selected = selected.genes,
                width='100%')
})

# multi select for Stimuli, reference table view
output$refTableStimuliSelectControl <- renderUI({
    stimuli.names <- setdiff(names(df.ref.table.qvalue),"variable")
    selectizeInput("ref_table_stimuli",
                   label = "Select Stimuli of interest",
                   choices = stimuli.names, 
                   multiple=TRUE,
                   selected = c('Null','FSL','pIC','LPS','FLA','GARD','ODN'),
                   options = list(placeholder = 'select a stimulus'))
})

#
# Venn Diagram select control
#
output$vennDiagramSelectControl <- renderUI({
    
    choices <- c("IFN-B / IFN-G / IL1B / TNFA" = "I3T",
                 "Cytokines / TLRs / Microbes" = "CTM")
    
    selectInput("venn_id", "Venn Diagram Setting",
                choices,
                selected=c('CTM'))
})

#
# mRNA signature clustering select control on stimuli
#
output$stimuliSelectionSignatures <- renderUI({
    names <- levels(df.data.long$StimulusName)
    selectizeInput("mrna_signature_stimuli",
                   label = "Select Stimuli for mRNA Signature:",
                   choices = names, 
                   multiple=TRUE,
                   selected = c('HKEC','HKSA','HKHP','HKLR','HKCA','BCG','IAV','SeV'),
                   options = list(placeholder = 'select a stimulus'))
}) 

#
# Stimuli Correlation
#
output$stimuliCorrelationmRNASelectControl <- renderUI({
    #names <- levels(df.data.long$mRNA)
    names <- unique(df.data.long$mRNA)
    selectizeInput("mrna_stimuli_correlations",
                   label = "Select mRNA",
                   choices = names, 
                   multiple=FALSE,
                   selected = "IRAK3")
})

#
# mRNA signature clustering select control on mRNA
#
output$mRNAsSelectionSignatures <- renderUI({
    names <- levels(df.data.long$mRNA)
    selectizeInput("mrna_signature_mRNAs",
                   label = "Select mRNAs to highlight:",
                   choices = names, 
                   multiple=TRUE,
                   selected = NULL,
                   options = list(placeholder = 'select mRNAs to highlight'))
}) 

#
# Metric MDS selection
#
output$metricSelectionSignatures <- renderUI({
    choices <- c("Max Normalization" = "MaxNorm",
                 "Z-score normalization" = "ZScoreNorm")
    
    selectInput("metricNorm", "Normalization Method",
                choices,
                selected=c('ZScoreNorm'))
})

output$fig1AButton <- renderUI({
    actionButton("pca.fig1A", "Fig.1A")
})

output$fig1BButton <- renderUI({
    actionButton("pca.fig1B", "Fig.1B")
})

output$figS4AButton <- renderUI({
    actionButton("pca.figS4A", "Fig.S4A")
})

output$figS4BButton <- renderUI({
    actionButton("pca.figS4B", "Fig.S4B")
})

output$figS5Button <- renderUI({
    actionButton("pca.figS5", "Fig.S5")
})

# boxplot
output$figSXXXButton <- renderUI({
    actionButton("boxplot.figSXXX", "Fig.SXXX")
})

output$fig2AButton <- renderUI({
    actionButton("boxplot.fig2A", "Fig.2A")
})

output$fig2CButton <- renderUI({
    actionButton("boxplot.fig2C", "Fig.2C")
})

output$fig3BButton <- renderUI({
    actionButton("boxplot.fig3B", "Fig.3B")
})