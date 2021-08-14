# Author details
## Script name: Figure 4.R
## Purpose of the script: To visualize the boxplot from the DNA methylation deconvolution data mentioned in the Das et al (2021), Epigenetics
## Author(s): Jyotirmoy Das, Ph.D.
## Date Created: 2019-03-21
## Date Last Modified: 2021-06-20
## Copyright statement: MIT open license
## Contact information: jyotirmoy.das@liu.se
## Please cite: 
## @article{das2021dna,
##  title={DNA methylome-based validation of induced sputum as an effective protocol to study lung immunity: construction of a classifier of pulmonary cell types},
##  author={Das, Jyotirmoy and Idh, Nina and Paues, Jacob and Sikkeland, Liv Ingunn Bjoner and Lerm, Maria},
##  journal={bioRxiv},
##  year={2021},
##  publisher={Cold Spring Harbor Laboratory}}
## Notes: The user needs to have the IDAT files which can be available upon request to the corresponding author.
## Please read the README file for the sample preparation of the analysis.

deconvolutionPlot <- function(deconvolutionPlot){
    #setwd("~/Desktop/Link to BackUp_All/BALdata/BALdata/Part2/DataSheets")
    uniqCpGs <- read.table("cell-specific-CpGs_ourStudy.txt", stringsAsFactors = F, header = T)
    BAL_HLADR <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors = F, header = T)
    BAL_CD3 <- read.table("myNorm_BAL_CD3.txt", stringsAsFactors = F, header = T)
    uniqCpGsHLADR <- uniqCpGs[1:594,]
    uniqCpGsCD3 <- uniqCpGs[595:2886,]

    suppressMessages(suppressWarnings(library(EpiDISH)))

    referenceData <- read.table("refDataB1_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData <- read.table("Stockholm_normalizedValue_GSE133062.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData.m <- as.matrix(methData)
    refData.m <- as.matrix(referenceData)
    out.lB1 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    # create individual plot, Sample P1
    boxplot(out.lB1, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB2_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    refData.m <- as.matrix(referenceData)
    out.lB2 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    # create individual plot, Sample P2
    boxplot(out.lB2, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB3_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    refData.m <- as.matrix(referenceData)
    out.lB3 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    # create individual plot, Sample P3
    boxplot(out.lB3, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB4_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    refData.m <- as.matrix(referenceData)
    out.lB4 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    # create individual plot, Sample P4
    boxplot(out.lB4, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB5_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    refData.m <- as.matrix(referenceData)
    out.lB5 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    # create individual plot, Sample P5
    boxplot(out.lB5, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataAVG_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    refData.m <- as.matrix(referenceData)
    out.lBAVG <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    # create individual plot, Sample avg
    boxplot(out.lBAVG, col = c("#FDD4D4", "#99E5E0"))

    # Combine all plots 
    outer = FALSE
    line = 1
    cex = 1.3
    adj  = 0.025
    par(mfrow=c(3,2), tcl=-0.5, family="Times New Roman", omi=c(0.3,0.3,0,0), 
        mar=c(2.3, 2.3, 1.8, 1.8), oma = c(2, 2, 1, 1))
    boxplot(out.lB1,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B1",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.lB2,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B2",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.lB3,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B3",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.lB4,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B4",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.lB5,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B5",cex.main=cex,col="black",font=1,line=line)
    #boxplot(out.lBAVG, col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    #title(outer=outer,adj=adj,main="Average",cex.main=cex,col="black",font=1,line=line)
    mtext("Cell proportion", side = 2, outer = T, at = 0.5, cex = 1.3)
    mtext("Cell types", side=1, outer=T, at=0.5, cex = 1.3)
}

deconvolutionPlot2 <- function(deconvolutionPlot2){
    #Figure 4b or Figure 5; Next part
    # Plan is to find how similar (in case of cell proportion) BAL and SI samples when BAL samples (CpGs) are used as
    # the reference dataset and SI samples (CpGs) are methData set

    referenceData <- read.table("refDataB1_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData <- read.table("TestDataB1_SI_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData.m <- as.matrix(methData)
    refData.m <- as.matrix(referenceData)
    out.testB1 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.testB1, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB2_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData <- read.table("TestDataB2_SI_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData.m <- as.matrix(methData)
    refData.m <- as.matrix(referenceData)
    out.testB2 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.testB2, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB3_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData <- read.table("TestDataB3_SI_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData.m <- as.matrix(methData)
    refData.m <- as.matrix(referenceData)
    out.testB3 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.testB3, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB4_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData <- read.table("TestDataB4_SI_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData.m <- as.matrix(methData)
    refData.m <- as.matrix(referenceData)
    out.testB4 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.testB4, col = c("#FDD4D4", "#99E5E0"))

    referenceData <- read.table("refDataB5_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData <- read.table("TestDataB5_SI_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
    methData.m <- as.matrix(methData)
    refData.m <- as.matrix(referenceData)
    out.testB5 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.testB5, col = c("#FDD4D4", "#99E5E0"))

    outer = FALSE
    line = 1
    cex = 1.3
    adj  = 0.025
    par(mfrow=c(3,2), tcl=-0.5, family="Times New Roman", omi=c(0.3,0.3,0,0), 
        mar=c(2.3, 2.3, 1.8, 1.8), oma = c(2, 2, 1, 1))
    boxplot(out.testB1,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B1",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.testB2,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B2",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.testB3,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B3",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.testB4,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B4",cex.main=cex,col="black",font=1,line=line)
    boxplot(out.testB5,  col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
    title(outer=outer,adj=adj,main="B5",cex.main=cex,col="black",font=1,line=line)
    mtext("Cell proportion", side = 2, outer = T, at = 0.5, cex = 1.3)
    mtext("Cell types", side=1, outer=T, at=0.5, cex = 1.3)
}

# Generating Correlation plot, Box plot and Violin plot
CellDensityBoxViolin <- function(CellDensityBoxViolin){
    suppressMessages(suppressWarnings(library(corrplot)))
    mB <- cor(data)
    mB1 <- as.matrix(B1_data)
    # Correlation plots
    corrplot(mB, method = "pie")
    corrplot(mB, method = "number")
    boxplot(mB)
    SI <- grepl("_SI_", colnames(data))
    colo1 <- ifelse(SI, "#FDD4D4", "#99E5E0")
    par(mar = c(6.1, 4.1, 4.1, 4.1), # change the margins
        lwd = 2, # increase the line thickness
        cex.axis = 1.2 # increase default axis label size
    )
    # Box plot 
    boxplot(mB, xaxt = "n", yaxt = "n", col = colo1, notch = T)
    axis(side = 1, labels = F)
    axis(side = 2)
    text(x = 1:length(data),
         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
         labels = colnames(data),
         col = colo2,
         family = "Times New Roman",
         font.size = 12,
         xpd = TRUE,
         adj = c(0.6, -3.5),
         offset = 1,
         ## Rotate the labels by 35 degrees.
         srt = 40,
         cex = 1.1)
    abline(h=0.9, col = "red", lty = 3)

    HLADR <- grepl("_HLADR", colnames(data))
    colo2 <- ifelse(HLADR, "#F97806", "#0BF488")
    par(mar = c(6.1, 4.1, 4.1, 4.1), # change the margins
        lwd = 2, # increase the line thickness
        cex.axis = 1.2 # increase default axis label size
    )
    # Violin plot
    vioplot(mB, col = colo1, xaxt = "n")
    text(x = 1:length(data),
         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
         labels = colnames(data),
         col = colo2,
         family = "Times New Roman",
         xpd = TRUE,
         adj = c(0.6, -3.5),
         offset = 1,
         ## Rotate the labels by 35 degrees.
         srt = 40,
         cex = 1.1)
    abline(h=0.9, col = "red", lty = 3)
}