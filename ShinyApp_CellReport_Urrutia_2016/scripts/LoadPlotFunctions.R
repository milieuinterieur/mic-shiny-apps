#
# LoadPlotFunctions.R
#
#

#
# getVennDiagram from t-test results
#

require(dplyr)
require(ggplot2)


#
# Facetted mRNA boxplots
#
#

facettedmRNABoxplots <- function(df.data = NULL, df.stimulus.color = NULL, compactPlotFlag = FALSE){
    
    cbPalette <- df.stimulus.color[match(as.character(levels(df.data$StimulusName)), df.stimulus.color[,1]),2]
    
    # Boxplot in ggplot2
    p <- ggplot(data=df.data, aes_string(x='StimulusName', y='value', colour='StimulusName'))
    p <- p + geom_boxplot(outlier.shape = NA)
    p <- p + geom_jitter()
    p <- p + facet_wrap(~ mRNA, ncol = 1, scales = ifelse(!compactPlotFlag,"free","free_y"))
    p <- p + scale_y_log10() + expand_limits(y=10)
    p <- p + scale_colour_manual(values = as.character(cbPalette))
    p <- p + theme_bw()
    p <- p + ylab("nCounts (log10)") + xlab("Stimulus Name")
    p <- p + theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size=16),strip.text.x = element_text(size=16))
    
    return(p)
}


#
# Facetted PCA Plot
#
facettedProjectedPcaPlot <- function(df.data=NULL, xName=NULL, yName=NULL, ellipse.stimuli=NULL, points.stimuli=NULL, stimuli.colors=NULL, 
                                     labels=TRUE, ellipse.linetype=1, xy.limits=c(-20,20)){
    
    require(ellipse)
    df.data.labels <- df.data %>% 
        filter(StimulusName %in% ellipse.stimuli) %>%
        group_by(StimulusName) %>%
        summarise_each(funs(mean))
    #summarise(meanPC1=mean(PC1), meanPC2=mean(PC2))
    
    names(df.data.labels)[1] <- "StimulusName2"
    
    #calculating ellipses
    require(ellipse)
    df_ell <- data.frame()
    df_points <- df.data
    names(df_points)[3:4] <- c('X','Y')
    for(g in ellipse.stimuli){
        df_ell <- rbind(df_ell, cbind(as.data.frame(with(df_points[df_points$StimulusName==g,], ellipse:::ellipse(cor(X, Y), 
                                                                                                                  scale=c(sd(X),sd(Y)), 
                                                                                                                  centre=c(mean(X),mean(Y)), level=0.9))),group=g))
    }
    if (length(points.stimuli)==0){
        p <- ggplot(data=NULL, aes_string(x=xName, y=yName, colour='StimulusName'))
        p <- p + theme_bw()
        p <- p + geom_vline(xintercept = 0, size=1) + geom_hline(yintercept=0,size=1)
        p <- p + geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=rep(1, times=length(ellipse.linetype)), linetype=ellipse.linetype)
        p <- p + scale_colour_manual(values = stimuli.colors)
        if (labels){
            p <- p + geom_text(data=df.data.labels, aes_string(x=xName, y=yName,label="StimulusName2", colour="StimulusName2", size=18))
        } else {
            p <- p + geom_text(data=subset(df.data.labels, StimulusName2 %in% c('TNFa','IL-1B','IFN-B','IFN-G')), aes_string(x=xName, y=yName,label="StimulusName2", colour="StimulusName2", size=18))
        }
        p <- p + theme(legend.position="none")
        p <- p + theme(panel.border = element_rect(linetype = 1, colour = "black", size=2))
        p <- p + scale_x_reverse() + scale_y_reverse()
        p <- p + coord_fixed(ratio = 1)#,xlim = xy.limits, ylim=xy.limits)
        p <- p + coord_equal()
    } else {
        p <- ggplot(data=df.data[df.data$StimulusName %in% points.stimuli,], aes_string(x=xName, y=yName, colour='StimulusName'))
        p <- p + geom_point(size=2.5)
        p <- p + theme_bw()
        p <- p + geom_vline(xintercept = 0, size=1) + geom_hline(yintercept=0,size=1)
        p <- p + geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1)
        p <- p + scale_colour_manual(values = stimuli.colors)
        p <- p + geom_text(data=df.data.labels, aes_string(x=xName, y=yName,label="StimulusName2", colour="StimulusName2", size=14))
        p <- p + facet_wrap(~ StimulusName, ncol = 1)
        p <- p + theme(legend.position="none")
        p <- p + theme(panel.border = element_rect(linetype = 1, colour = "black", size=2))
        p <- p + theme(strip.text.x = element_text(colour = "black", face='bold', angle = 0, size = 14))
        p <- p + scale_x_reverse() + scale_y_reverse() 
    }
    
    return(p)
}



getVennDiagramFromTtest <- function(df.ttest.data=NULL, stimuli=NULL, colors=NULL, alpha.coef=0.5){
    require(VennDiagram)
    
    df.ttest.data <- as.data.frame(df.ttest.data)
    
    #get stimulation names
    names.stimulation <- stimuli
    
    #convert data.frame into list
    list.diff.expression <- list()
    
    # extract expressed gene from each stimulation
    for(stim in names.stimulation){
        v.tmp <- as.character(df.ttest.data[df.ttest.data$StimulusName.y==stim & df.ttest.data$adj.p.value<0.01,2,drop=TRUE])
        #print(v.tmp)
        list.diff.expression[[paste(stim," \n(",length(v.tmp),")",sep="")]] <- v.tmp
        
    }
    
    # build Venn Diagram Grid
    venn.plot.grid <- venn.diagram(
        x = list.diff.expression,
        filename = NULL,
        main = "", #paste("Total number of genes under consideration: ",length(unique(df.ttest.data[,2])),sep=""), #,
        cat.col = colors,
        fill = colors,
        fontface = 15,
        alpha = 0.7,
          main.col='black',
        label.col = "white",
        #alpha = alpha.coef,
        cex = 2.5,
        cat.cex = 1.8,
        cat.dist = c(0.25,0.25,0.2,0.2)
        #cat.dist = rep(x = 0.05, times = length(names.stimulation))
    );
    #grid.draw(venn.plot);
    return(venn.plot.grid)
}

getVennDiagramFromDataFrame <- function(df.data=NULL, colors=NULL, alpha.coef=0.5){
  require(VennDiagram)
  
  df.data <- read.table(file='data//Cytokines_TLRs_Bugs_geneLists.txt', header=TRUE, sep='\t', stringsAsFactors = FALSE)
  #get stimulation names
  names.groups <- names(df.data)
  alpha.coef <- 0.5
  #convert data.frame into list
  list.diff.expression <- as.list(df.data)
  
  colors <- c('red','blue','green')
  # build Venn Diagram Grid
  venn.plot.grid <- venn.diagram(
    x = list.diff.expression,
    filename = NULL,
    main = "Main Title", #,
    #cat.col = colors,
    fill = colors,
    alpha = alpha.coef
    #cex = 2.5,
    #cat.cex = 2.5,
    #cat.pos = 0,
    #cat.dist = rep(x = 0.05, times = length(names.stimulation))
  );
  #grid.draw(venn.plot);
  return(venn.plot.grid)
}


getMDSPlot <- function(df.data=NULL, selected.stimuli=NULL, selected.mRNAs=NULL){
    
    #for(stimu in selected.stimuli){
        df.curr <- subset(df.data, StimulusName %in% selected.stimuli)
        cor.mRNAs <- cor(df.curr[,selected.mRNAs], 
                         method="spearman",use="pairwise.complete.obs")
        cor.mRNAs.diss <- 1-abs(cor.mRNAs)
        
        fit <- cmdscale(cor.mRNAs.diss,eig=TRUE, k=2) 
        
        # plot solution
        df.mds.data <- data.frame(xPoint=fit$points[,1], 
                                  yPoint=fit$points[,2], 
                                  mRNAs=row.names(fit$points))
        
        p <- ggplot(data=df.mds.data, aes_string(x='xPoint', y='yPoint', label='mRNAs'))
        p <- geom_text(size=10)
        p <- p + theme_bw()
        return(p)
    #}
}

#
# get CorrelationMatrix
#
plotGgplot2CorrMatrix <- function(df.data=NULL){#, selected.stimuli=NULL){
  
  require(reshape2)
  
  df.x <- df.data #cbind(names=rownames(df.data), df.data)
  molten.rbm <- melt(df.x) # cbind(melt(c), stars) 
  names(molten.rbm) <- c("M1", "M2", "corr") #, "pvalue")
  #molten.rbm$M1 <- (as.character(molten.rbm$M1))
  #molten.rbm$M2 <- (as.character(molten.rbm$M2))
  
  #define each triangle of the plot matrix and the diagonal (mi.ids)
  rbm.ids <- subset(molten.rbm, M1 == M2)
  rbm.lower <- subset(molten.rbm[lower.tri(df.x),], M1 != M2)
  rbm.upper <- subset(molten.rbm[upper.tri(df.x),], M1 != M2)
  
  if (!all(is.na(rbm.lower$corr))){
    p <- ggplot(molten.rbm, aes(x=M1, M2, fill=corr))# + theme_bw() # + geom_tile()    
    # now plot just these values, adding labels (geom_text) for the names and th values
    p <- p + geom_tile(data=rbm.lower)#, width=0.8, height=0.8)
    p <- p + geom_text(data=rbm.lower, size=6, aes(label=round(corr,2)))#paste(round(corr,3), pvalue)))
    p <- p + geom_text(data=rbm.ids, size=6, aes(label=M2, colour="black", angle=0, hjust= 0.7)) 
    
    meas <- as.character(rev(unique(molten.rbm$M2)))
    p <- p + scale_colour_identity()  
    p <- p + scale_fill_gradientn(colours= c("red", "#F0F0F0","#F0F0F0", "#F0F0F0","blue"), limits =c(1,-1))
    
    p <- p + scale_x_discrete(limits=meas[length(meas):1]) 
    p <- p + scale_y_discrete(limits=meas) 
    
    p <- p + xlab(NULL) + ylab(NULL) 
    
    p <- p + theme_bw()
    p <- p + theme(legend.position = "none")
    p <- p + theme(axis.title=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
    p <- p + theme(plot.title = element_text(size = rel(2)))
    
    return(p)
  } else {
    print(paste("Empty Correlation Matrix ") )
    return(NULL)
  }
}