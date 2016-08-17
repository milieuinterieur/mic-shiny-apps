#
# LoadDataFunctions.R
# v0.4
#


#
# Load Stimuli Info
#

getStimulusList <- function(){
    df.stimuli.selected <- read.table(file='data//TC_stimuli_table_v6.txt', header=TRUE, sep='\t')
    return(df.stimuli.selected)
}

#
# Load Stimuli Color Palette
#
getStimulusColorPalette <- function(){
    df.stimuli.selected <- read.table(file='data//TC_stimuli_table_v6.txt', header=TRUE, sep='\t')
    return(df.stimuli.selected[, c('StimulusName','Color')])
}

#
# Optimized gene list on I3T cytokines
# following t-test + projection score optimisation
# source Qlucore
getI3Tgenes <- function(){
    require(stringr)
    df.genes <- read.table(file="data//Nanostring25_SelectedNanostringProbesList_q0-001_ps0-62_192genes_17062015.txt", header=TRUE, sep='\t',stringsAsFactors = FALSE)
    genes <- as.character(df.genes$selected.probes)
    genes <- str_replace(genes, pattern = "\\/", replacement = ".")
    genes <- str_replace(genes, pattern = "-", replacement = ".")
    return(genes)
}

get202genes <- function(){
    df.genes.202 <- read.table(file="data/I3T_202_genes_12062015.txt", header=TRUE, sep='\t',stringsAsFactors = FALSE)
    genes <- as.character(df.genes.202$I3T_202_genes)
    genes <- str_replace(genes, pattern = "\\/", replacement = ".")
    genes <- str_replace(genes, pattern = "-", replacement = ".")
    return(genes)
}

# Top 40 genes from I3T selection
getTop40Genes <- function(){   
    ifnb.top10.genes <- c('CCL8','LAG3','MX1','IFIH1','BST2','IFIT2','IFITM1','IRF7','CCR1','STAT2')
    il1b.top10.genes <- c('IL6','CXCL2','CXCL1','MARCO','CCL7','IL1B','IL1A','CXCL13','CD163','CCL20')
    ifng.top10.genes <- c('CXCL9','HLA.DPA1','TLR8','C1QB','RARRES3','HLA.DMA','PDCD1LG2','FCGR1A.B','SOCS1','IRF1')
    tnfa.top10.genes <- c('CD83','NFKB2','TNFAIP3', 'C3','IL18','CCL23','EBI3','CLEC6A','TNFRSF4','CLEC4E')
    
    top40.genes <- c(ifnb.top10.genes, il1b.top10.genes, ifng.top10.genes, tnfa.top10.genes)
    
    return(top40.genes)  
}

#
# Optimized gene list on I3T cytokines
# following t-test + projection score optimisation
# source Qlucore
getSelectedProbesgenes <- function(){
    require(stringr)
    df.genes <- read.table(file="data/Nanostring25_SelectedNanostringProbesList_572genes_17062015.txt", header=TRUE, sep='\t',stringsAsFactors = FALSE)
    genesQC <- as.character(df.genes$selected.probes)
    genesQC <- str_replace(genesQC, pattern = "\\/", replacement = ".")
    genesQC <- str_replace(genesQC, pattern = "-", replacement = ".")
    return(genesQC)
}

# get Adjusted q values of Paired T-Test
getQValues <- function(){
    
    df.qvalues <- read.table(file='data/Nanostring25_ReferenceValues_Stats_qValue_2015-08-15.tsv', 
                             stringsAsFactors = FALSE,
                             sep='\t',
                             header=TRUE)
    return(df.qvalues)
}

# get Fold Change Expression
getFCValues <- function(){
    
    df.fcvalues <- read.table(file='data/Nanostring25_ReferenceValues_Stats_FC_2015-08-15.tsv', 
                              stringsAsFactors = FALSE,
                              sep='\t',
                              header=TRUE)
    return(df.fcvalues)
}

# get Median Expression
getMedianValues <- function(){
    
    df.fcvalues <- read.table(file='data/Nanostring25_ReferenceValues_Stats_Median_2015-08-15.tsv', 
                              stringsAsFactors = FALSE,
                              sep='\t',
                              header=TRUE)
    return(df.fcvalues)
}

# get CV Expression
getCvValues <- function(){
    
    df.fcvalues <- read.table(file='data/Nanostring25_ReferenceValues_Stats_CV_2015-08-15.tsv', 
                              stringsAsFactors = FALSE,
                              sep='\t',
                              header=TRUE)
    return(df.fcvalues)
}


getNanoStringFinalDataResubmission <- function(selectedProbes=TRUE){
    
    require(stringr)
    #require(raster)
    require(dplyr)
    require(reshape2)
    require(ggplot2)
    
    # Load All Stimulus Coding information (40 stimulations)
    df.stimuli.info <- read.table(file='data/TC_stimuli_table_v6.txt', 
                                  stringsAsFactors = FALSE,
                                  sep='\t',
                                  header=TRUE)
    # keep same stimuli order as StimulusId
    df.stimuli.info$StimulusName <- factor(df.stimuli.info$StimulusName,
                                           levels = df.stimuli.info$StimulusName)
    
    # Load NanoString data from TSV Export
    #df.nano <- read.table(file='data/20150527_1000RNA_25Donors_HKNormalized_SampleRow.csv',
    df.nano <- read.table(file='data/20150929_1000RNA_25Donors_HKNormalized_SampleRow.csv',
                          stringsAsFactors = FALSE,
                          sep=';',
                          header=TRUE)
    #v.all.mRNAs <- names(df.nano)[c(17:610)]
    v.all.mRNAs <- names(df.nano)[c(17:599)]
    
    # Normalization Genes: RPL19, TBP, POLR2A, HPRT1
    df.nano <- df.nano[1:1000,c('Sample_ID',v.all.mRNAs)]
     
    #Generate SampleID, DonorID, StimulationID Columns
    #
    df.nano <- df.nano %>%
                mutate(Donor.ID=str_sub(Sample_ID, start = 2, end = 5), 
                       StimulationId=as.numeric(str_sub(Sample_ID, start = 8, end = 9)))
    
    # Essential meta info columns for further analysis + RNA columns
    v.meta.info <- c('Sample_ID','Donor.ID','StimulationId')
    #v.mRNA.names <- getNanoStringGeneListV2()
    if (selectedProbes){
      v.mRNA.names <- getSelectedProbesgenes()
    } else {
      v.mRNA.names <- v.all.mRNAs
    }
    
    # convert Donor Ids into factors
    df.nano$Donor.ID <- factor(df.nano$Donor.ID)
    # merge both tables to get StimulusName
    df.nano <- merge(df.stimuli.info[,c('StimulusId','StimulusName', 'Color')], 
                     df.nano[,c(v.meta.info,v.mRNA.names)],  
                     by.x=c('StimulusId'), by.y=c('StimulationId'), all.y=TRUE)
    
    return(df.nano)
}

#
# Optimized 44 gene list on I3T cytokines
# source Qlucore
get44genes <- function(){
  require(stringr)
  df.genes <- read.table(file="data/Union of all 4 signatures.txt", header=FALSE, sep='\t',stringsAsFactors = FALSE)
  genes <- as.character(df.genes$V1)
  genes <- str_replace(genes, pattern = "\\/", replacement = ".")
  genes <- str_replace(genes, pattern = "-", replacement = ".")
  return(genes)
}

#
# Get correlations
#

getStimuliCorrelations <- function(){
    df.stimuli.corr <- read.table(file='data/Nanostring25Resubmission_SpearmanCorrelations_20151025.tsv', sep='\t', header=TRUE)
    return(df.stimuli.corr)
}

#
# Get paired t-test results
# 
getTTestsResults <- function(){
    df.ttest.results <- read.table(file='data/Nanostring25Resubmission_ttests_20151025.tsv', sep='\t',header=TRUE)
    return(df.ttest.results)
}