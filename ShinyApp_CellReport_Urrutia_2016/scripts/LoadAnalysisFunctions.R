#
#
# Load Analysis Functions
#
#

#
# Helper function to compute residuals following a linear regression
#
regressOnFormula<-function(n, data, formula=NULL) { 
    currFormula <- as.formula(paste(n, " ~ ", formula, sep=""))
    #print(currFormula)
    residus <- resid(lm(currFormula, data=data, na.action=na.exclude))
    return(residus)
}

#
# Returns residuals as a data.frame following a linear regression
#
getRegressedOutDF <- function(df.data=NULL, regression.formula=NULL, cols.to.regress=NULL, cols.to.keep=NULL){
    
    # get residus from all regressions
    l.data.regressed <- lapply(cols.to.regress, regressOnFormula, data = df.data, formula=regression.formula)
    # aggregate into data.frame
    df.data.regressed <- as.data.frame(do.call(cbind, l.data.regressed))
    names(df.data.regressed) <- cols.to.regress
    # bind with cols.to.keep from original data
    df.data.regressed <- cbind(df.data[,cols.to.keep], df.data.regressed)
    
    return(df.data.regressed)
}

#
# Compute PCA and Projections
#

getCoordinatesFromPcaAndProjections <- function(nbDim=2, df.data=NULL, pca.stimuli=NULL, projected.stimuli=NULL){
    
    require(FactoMineR)
    
    # Build composite data.frame
    # sort out data.frame (selected individual and non ones)
    df.data.selected.PCA <- df.data %>%
        filter(StimulusName %in% pca.stimuli) 
    
    # compute PCA
    if(!is.null(projected.stimuli)){
        # select individuals for PCA projection
        df.data.projected.PCA <- df.data %>% filter(StimulusName %in% projected.stimuli)
        # bind them
        df.data.PCA <- rbind(df.data.selected.PCA,df.data.projected.PCA)
        first.ind.supp <- dim(df.data.selected.PCA)[1] + 1
        last.ind.supp <- dim(df.data.projected.PCA)[1] + first.ind.supp -1
        
    pca.res <- PCA(X = df.data.PCA, 
                   ncp = nbDim, scale.unit = TRUE, quali.sup = 1:4, 
                   ind.sup = first.ind.supp:last.ind.supp, graph = FALSE)
    
    pca.coord <- pca.res$ind$coord
    pca.coord.indsupp <- pca.res$ind.sup$coord
    df.scatterplot <- df.data.PCA[,c('StimulusName','Color','Donor.ID')]
    for(i in 1:nbDim){
        df.scatterplot <- cbind(df.scatterplot, c(pca.coord[,i],pca.coord.indsupp[,i]))
    }
    } else {
        df.data.PCA <- df.data.selected.PCA
        
        pca.res <- PCA(X = df.data.PCA, 
                       ncp = nbDim, scale.unit = TRUE, quali.sup = 1:4, 
                       ind.sup = NULL, graph = FALSE)
        
        pca.coord <- pca.res$ind$coord
        df.scatterplot <- df.data.PCA[,c('StimulusName','Color','Donor.ID')]
        for(i in 1:nbDim){
            df.scatterplot <- cbind(df.scatterplot, pca.coord[,i])
        }  
    }
    
    
    names(df.scatterplot) <- c('StimulusName', 'Color', 'Donor.ID', paste("PC",1:nbDim,sep=""))
    
    return(df.scatterplot)
}

getPcaAndProjections <- function(nbDim=2, df.data=NULL, pca.stimuli=NULL, projected.stimuli=NULL){
    
    require(FactoMineR)
    
    # Build composite data.frame
    # sort out data.frame (selected individual and non ones)
    df.data.selected.PCA <- df.data %>%
        filter(StimulusName %in% pca.stimuli) 
    # select individuals for PCA projection
    df.data.projected.PCA <- df.data %>%
        filter(StimulusName %in% projected.stimuli)
    # bind them
    df.data.PCA <- rbind(df.data.selected.PCA,df.data.projected.PCA)
    first.ind.supp <- dim(df.data.selected.PCA)[1] + 1
    last.ind.supp <- dim(df.data.projected.PCA)[1] + first.ind.supp -1
    
    # select on genes ????
    pca.res <- PCA(X = df.data.PCA, 
                   ncp = nbDim, scale.unit = TRUE, quali.sup = 1:4, 
                   ind.sup = first.ind.supp:last.ind.supp, graph = FALSE)
    return(pca.res)
}

#
# MDS ready dataset to compare stimulus signatures
# log10 centered data
#
getMDSReadyData <- function(df.data.long=NULL){
    
    # Compute median stats
    df.data.long.median.stats <- df.data.long %>%
        group_by(StimulusName, mRNA) %>%
        summarise(median = median(value, na.rm = TRUE)) %>%
        spread(StimulusName, median)
    
    # apply log10 and substract Null stimulus
    df.data.long.median.stats.log10.centered <- as.data.frame(apply((log10(df.data.long.median.stats[,-1])), 2, function(y) y - (log10(df.data.long.median.stats[,2]))))
    names(df.data.long.median.stats.log10.centered) <- names(df.data.long.median.stats[,-1])
    #head(df.nano.filter.melted.median.stats.centered)
    row.names(df.data.long.median.stats.log10.centered) <- df.data.long.median.stats$mRNA
    
    return(df.data.long.median.stats.log10.centered)
}
