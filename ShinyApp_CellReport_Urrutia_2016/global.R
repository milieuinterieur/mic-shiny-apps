#
# Shiny Nanostring25
# Global.R
#

library(shiny)
library(BH)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(RColorBrewer)
library(markdown)
library(stringr)
library(reshape2)
library(corrplot)

#Suppress Warnings
options(warn=-1)

source('scripts/LoadDataFunctions.R')
source('scripts/LoadAnalysisFunctions.R')
source('scripts/LoadPlotFunctions.R')

df.stimulus.color <- getStimulusColorPalette()
I3T.stimuli <- c('IFN-B','IFN-G','IL-1B','TNFa')
TLR.stimuli <- c('FSL','pIC','LPS','FLA','GARD','r848','ODN')
BUG.stimuli <- c('HKHP','HKSA','HKLR','HKEC','BCG','HKCA','IAV','SeV')
selected.stimuli <- c("Null", I3T.stimuli, TLR.stimuli, BUG.stimuli,
                      "WGP", "LAM", "CPPD", "IFN-A") #,"DAP")
# Get Data
I3T.genes <- getI3Tgenes()
I3T.202.genes <- get202genes()
top40.genes <- getTop40Genes()


df.data.wide.resubmission <- getNanoStringFinalDataResubmission() %>% 
                                filter(StimulusName %in% selected.stimuli)
df.data.wide.resubmission$StimulusName <- factor(as.character(df.data.wide.resubmission$StimulusName),
                                                 levels = selected.stimuli)
selected.mRNAs <- names(df.data.wide.resubmission)[6:577]
df.data.long <- df.data.wide.resubmission %>% gather("mRNA","value",
                                                     one_of(selected.mRNAs))
df.data.wide.resubmission[,selected.mRNAs] <- log2(df.data.wide.resubmission[,selected.mRNAs])
df.data.wide <- df.data.wide.resubmission
top44.genes <- get44genes()
df.ref.table.qvalue <- getQValues()
df.stimuli.corr <- getStimuliCorrelations() %>% select(-SEB)
df.ttests.results <- getTTestsResults() %>% filter(StimulusName %in% selected.stimuli)

# Donor Regression on 192 genes
df.data.wide.DonorRegressed <- getRegressedOutDF(df.data=df.data.wide, 
                                                 regression.formula="Donor.ID", 
                                                 cols.to.regress=selected.mRNAs, 
                                                 cols.to.keep=c('StimulusName','Color','Donor.ID'))
