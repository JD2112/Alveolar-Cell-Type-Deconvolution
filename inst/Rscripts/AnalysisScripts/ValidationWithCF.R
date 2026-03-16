# Author details
## Script name: StatiscalAnalysis.R
## Purpose of the script: Different statiscal calculations from the DNA methylation deconvolution data mentioned in the Das et al (2021), Epigenetics
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
### New analysis with Cystic Fibrosis data
### Compare the BAL samples with CF data to find out the location
### of the sample collection
### 1. Get the regression analysis
### 2. Calculate the residual values
### 3. Select those CpGs who are at distant from the regression line (here: >\0.3\)
### 4. Use the identified CpGs as the training data to locate the lung position

ValidationWithCysticFibrosis <- function(ValidationWithCysticFibrosis){
    #load("~/Desktop/Link to BackUp_All/BALdata/BALdata/ServerData/HLADR/.RData")
    # collect residual information 
    b5_resd <- cbind(B5_CD3_1$data$ProbeI, mod1$residuals)
    colnames(b5_resd) <- c("ProbeID", "residuals")
    b5_resd <- cbind(B5_CD3_1$data, mod1$residuals)
    colnames(b5_resd)[12] <- "residuals"
    b5_resd_posMax <- b5_resd[b5_resd$residuals > 0.3,]
    b5_resd_negMax <- b5_resd[b5_resd$residuals < -0.3,]
    b5_resd_03 <- rbind(b5_resd_posMax, b5_resd_negMax)
    # Load Cystic Fibrosis dataset
    sampleTableCF <- read.table("SampleTable_GSE132547.txt", 
            stringsAsFactors = F, 
            header = F)
    colnames(sampleTableCF) <- c("ID", "POSITION")
    colnames(CFdata)[2:25] <- sampleTableCF$POSITION
    newcols <- sapply(c("RLL$", "RUL$"), function(x) rowMeans(CFdata[grep(x, names(CFdata))]))
    newcols1 <- setNames(cbind(CFdata[1], newcols), c(names(CFdata)[1], "RLL", "RUL"))
    # Perform the analysis
    suppressMessages(suppressWarnings(library(EpiDISH)))
    referenceData <- newcols1
    names(referenceData) <- newcols1[,1]
    rownames(referenceData) <- referenceData[,1]
    referenceData <- referenceData[,2:3]
    refData.m <- as.matrix(referenceData)
    methData <- b5_resd_03
    rownames(methData) <- methData[,1]
    methData <- methData[, 2:11]
    methData.m <- as.matrix(methData)
    out.l5 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l1, col = c("#FDD4D4", "#99E5E0"))
    #### Analysis with HLA-DR file
    # Load the data file
    hladr <- read.table("../ServerData/HLADR/myNorm_BAL_HLADR.txt", 
                stringsAsFactors = F, 
                header = T)
    # Perform the linear regression analysis
    lm_B1 <- lm(B1_SPU_HLDADR ~ B1_LAV_HLDADR, data = hladr)
    lm_B2 <- lm(B2_SPU_HLDADR ~ B2_LAV_HLDADR, data = hladr)
    lm_B3 <- lm(B3_SPU_HLDADR ~ B3_LAV_HLDADR, data = hladr)
    lm_B4 <- lm(B4_SPU_HLDADR ~ B4_LAV_HLDADR, data = hladr)
    lm_B5 <- lm(B5_SPU_HLDADR ~ B5_LAV_HLDADR, data = hladr)
    # Choose the residual data and combined with the CpG identifiers
    b1_resd <- cbind(hladr$ProbeID, as.data.frame(lm_B1$residuals))
    # Change the column names of the residual file and select residual > |0.3| CpGs
    colnames(b1_resd) <- c("ProbeID", "residuals")
    colnames(b1_resd) <- "residuals"
    b1_resd_posMax <- b1_resd[b1_resd$residuals > 0.3,]
    b1_resd_negMax <- b1_resd[b1_resd$residuals < -0.3,]
    b1_resd_03 <- rbind(b1_resd_posMax, b1_resd_negMax)
    b1_resd$ProbeID <- rownames(b1_resd)
    rownames(b1_resd) <- NULL
    b1_resd <- b1_resd[c(2,1)]
    b1_resd_posMax <- b1_resd[b1_resd$residuals > 0.3,]
    b1_resd_negMax <- b1_resd[b1_resd$residuals < -0.3,]
    b1_resd_03 <- rbind(b1_resd_posMax, b1_resd_negMax)
    b1_resd_03_data <- merge(b1_resd_03, hladr, by = "ProbeID")
    b1_resd_03_data_1 <- b1_resd_03_data[c(1,3,4)]
    rownames(b1_resd_03_data_1) <- b1_resd_03_data_1[,1]
    b1_resd_03_data_1 <- b1_resd_03_data_1[,2:3]
    methDataB1 <- b1_resd_03_data_1
    methDataB1.m <- as.matrix(methDataB1)
    out.l1 <- epidish(beta.m = methDataB1.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l1, col = c("#FDD4D4", "#99E5E0"))
    out.l1
    # Sample P2 analysis
    b2_resd <- cbind(hladr$ProbeID, as.data.frame(lm_B2$residuals))
    colnames(b2_resd) <- c("ProbeID", "residuals")
    rownames(b2_resd) <- NULL
    b2_resd_posMax <- b2_resd[b2_resd$residuals > 0.3,]
    b2_resd_negMax <- b2_resd[b2_resd$residuals < -0.3,]
    b2_resd_03 <- rbind(b2_resd_posMax, b2_resd_negMax)
    b2_resd_03_data <- merge(b2_resd_03, hladr, by = "ProbeID")
    b2_resd_03_data_1 <- b2_resd_03_data[c(1,5,6)]
    rownames(b2_resd_03_data_1) <- b2_resd_03_data_1[,1]
    b2_resd_03_data_1 <- b2_resd_03_data_1[,2:3]
    methDataB2 <- b2_resd_03_data_1
    methDataB2.m <- as.matrix(methDataB2)
    out.l2 <- epidish(beta.m = methDataB2.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l2, col = c("#FDD4D4", "#99E5E0"))
    out.l2
    # Sample P3 analysis
    b3_resd <- cbind(hladr$ProbeID, as.data.frame(lm_B3$residuals))
    colnames(b3_resd) <- c("ProbeID", "residuals")
    rownames(b3_resd) <- NULL
    b3_resd_posMax <- b3_resd[b3_resd$residuals > 0.3,]
    b3_resd_negMax <- b3_resd[b3_resd$residuals < -0.3,]
    b3_resd_03 <- rbind(b3_resd_posMax, b3_resd_negMax)
    b3_resd_03_data <- merge(b3_resd_03, hladr, by = "ProbeID")
    b3_resd_03_data_1 <- b3_resd_03_data[c(1,7,8)]
    rownames(b3_resd_03_data_1) <- b3_resd_03_data_1[,1]
    b3_resd_03_data_1 <- b3_resd_03_data_1[,2:3]
    methDataB3 <- b3_resd_03_data_1
    methDataB3.m <- as.matrix(methDataB3)
    out.l3 <- epidish(beta.m = methDataB3.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l3, col = c("#FDD4D4", "#99E5E0"))
    out.l3
    # Sample P4 analysis
    b4_resd <- cbind(hladr$ProbeID, as.data.frame(lm_B4$residuals))
    colnames(b4_resd) <- c("ProbeID", "residuals")
    rownames(b4_resd) <- NULL
    b4_resd_posMax <- b4_resd[b4_resd$residuals > 0.3,]
    b4_resd_negMax <- b4_resd[b4_resd$residuals < -0.3,]
    b4_resd_03 <- rbind(b4_resd_posMax, b4_resd_negMax)
    b4_resd_03_data <- merge(b4_resd_03, hladr, by = "ProbeID")
    b4_resd_03_data_1 <- b4_resd_03_data[c(1,9,10)]
    rownames(b4_resd_03_data_1) <- b4_resd_03_data_1[,1]
    b4_resd_03_data_1 <- b4_resd_03_data_1[,2:3]
    methDataB4 <- b4_resd_03_data_1
    methDataB4.m <- as.matrix(methDataB4)
    out.l4 <- epidish(beta.m = methDataB4.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l4, col = c("#FDD4D4", "#99E5E0"))
    out.l4
    # Smaple P5 analysis
    b5_resd <- cbind(hladr$ProbeID, as.data.frame(lm_B5$residuals))
    colnames(b5_resd) <- c("ProbeID", "residuals")
    rownames(b5_resd) <- NULL
    b5_resd_posMax <- b5_resd[b5_resd$residuals > 0.3,]
    b5_resd_negMax <- b5_resd[b5_resd$residuals < -0.3,]
    b5_resd_03 <- rbind(b5_resd_posMax, b5_resd_negMax)
    b5_resd_03_data <- merge(b5_resd_03, hladr, by = "ProbeID")
    b5_resd_03_data_1 <- b5_resd_03_data[c(1,11,12)]
    rownames(b5_resd_03_data_1) <- b5_resd_03_data_1[,1]
    b5_resd_03_data_1 <- b5_resd_03_data_1[,2:3]
    methDataB5 <- b5_resd_03_data_1
    methDataB5.m <- as.matrix(methDataB5)
    out.l5 <- epidish(beta.m = methDataB5.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l5, col = c("#FDD4D4", "#99E5E0"))
    out.l5
    # CD3 data analysis
    cd3 <- read.table("../ServerData/CD3/myNorm_BAL_CD3.txt", 
                stringsAsFactors = F, 
                header = T)

    # Smaple P5 analysis
    b5_resd <- cbind(hladr$ProbeID, as.data.frame(lm_B5$residuals))
    colnames(b5_resd) <- c("ProbeID", "residuals")
    rownames(b5_resd) <- NULL
    b5_resd_posMax <- b5_resd[b5_resd$residuals > 0.3,]
    b5_resd_negMax <- b5_resd[b5_resd$residuals < -0.3,]
    b5_resd_03 <- rbind(b5_resd_posMax, b5_resd_negMax)
    b5_resd_03_data <- merge(b5_resd_03, hladr, by = "ProbeID")
    b5_resd_03_data_1 <- b5_resd_03_data[c(1,11,12)]
    rownames(b5_resd_03_data_1) <- b5_resd_03_data_1[,1]
    b5_resd_03_data_1 <- b5_resd_03_data_1[,2:3]
    methDataB5 <- b5_resd_03_data_1
    methDataB5.m <- as.matrix(methDataB5)
    out.l5 <- epidish(beta.m = methDataB5.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.l5, col = c("#FDD4D4", "#99E5E0"))
    out.l5
    # Perform the linear regression analysis
    lm_cd3B1 <- lm(B1_SPU_CD3 ~ B1_LAV_CD3, data = cd3)
    lm_cd3B2 <- lm(B2_SPU_CD3 ~ B2_LAV_CD3, data = cd3)
    lm_cd3B3 <- lm(B3_SPU_CD3 ~ B3_LAV_CD3, data = cd3)
    lm_cd3B4 <- lm(B4_SPU_CD3 ~ B4_LAV_CD3, data = cd3)
    lm_cd3B5 <- lm(B5_SPU_CD3 ~ B5_LAV_CD3, data = cd3)
    # Choose the residual data and combined with the CpG identifiers
    cd3b1_resd <- cbind(cd3$ProbeID, as.data.frame(lm_cd3B1$residuals))
    colnames(cd3b1_resd) <- c("ProbeID", "residuals")
    rownames(cd3b1_resd) <- NULL
    cd3b1_resd_posMax <- cd3b1_resd[cd3b1_resd$residuals > 0.3,]
    cd3b1_resd_negMax <- cd3b1_resd[cd3b1_resd$residuals < -0.3,]
    cd3b1_resd_03 <- rbind(cd3b1_resd_posMax, cd3b1_resd_negMax)
    cd3b1_resd_03_data <- merge(cd3b1_resd_03, cd3, by = "ProbeID")
    cd3b1_resd_03_data_1 <- cd3b1_resd_03_data[c(1,3,4)]
    rownames(cd3b1_resd_03_data_1) <- cd3b1_resd_03_data_1[,1]
    cd3b1_resd_03_data_1 <- cd3b1_resd_03_data_1[,2:3]
    methDatacd3b1 <- cd3b1_resd_03_data_1
    methDatacd3b1.m <- as.matrix(methDatacd3b1)
    out.lcd3b1 <- epidish(beta.m = methDatacd3b1.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.lcd3b1, col = c("#FDD4D4", "#99E5E0"))
    out.lcd3b1
    # Choose the residual data and combined with the CpG identifiers
    cd3b2_resd <- cbind(cd3$ProbeID, as.data.frame(lm_cd3B2$residuals))
    colnames(cd3b2_resd) <- c("ProbeID", "residuals")
    rownames(cd3b2_resd) <- NULL
    cd3b2_resd_posMax <- cd3b2_resd[cd3b2_resd$residuals > 0.3,]
    cd3b2_resd_negMax <- cd3b2_resd[cd3b2_resd$residuals < -0.3,]
    cd3b2_resd_03 <- rbind(cd3b2_resd_posMax, cd3b2_resd_negMax)
    cd3b2_resd_03_data <- merge(cd3b2_resd_03, cd3, by = "ProbeID")
    cd3b2_resd_03_data_1 <- cd3b2_resd_03_data[c(1,5,6)]
    rownames(cd3b2_resd_03_data_1) <- cd3b2_resd_03_data_1[,1]
    cd3b2_resd_03_data_1 <- cd3b2_resd_03_data_1[,2:3]
    methDatacd3b2 <- cd3b2_resd_03_data_1
    methDatacd3b2.m <- as.matrix(methDatacd3b2)
    out.lcd3b2 <- epidish(beta.m = methDatacd3b2.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.lcd3b2, col = c("#FDD4D4", "#99E5E0"))
    out.lcd3b2
    # Choose the residual data and combined with the CpG identifiers
    cd3b3_resd <- cbind(cd3$ProbeID, as.data.frame(lm_cd3B3$residuals))
    colnames(cd3b3_resd) <- c("ProbeID", "residuals")
    rownames(cd3b3_resd) <- NULL
    cd3b3_resd_posMax <- cd3b3_resd[cd3b3_resd$residuals > 0.3,]
    cd3b3_resd_negMax <- cd3b3_resd[cd3b3_resd$residuals < -0.3,]
    cd3b3_resd_03 <- rbind(cd3b3_resd_posMax, cd3b3_resd_negMax)
    cd3b3_resd_03_data <- merge(cd3b3_resd_03, cd3, by = "ProbeID")
    cd3b3_resd_03_data_1 <- cd3b3_resd_03_data[c(1,7,8)]
    rownames(cd3b3_resd_03_data_1) <- cd3b3_resd_03_data_1[,1]
    cd3b3_resd_03_data_1 <- cd3b3_resd_03_data_1[,2:3]
    methDatacd3b3 <- cd3b3_resd_03_data_1
    methDatacd3b3.m <- as.matrix(methDatacd3b2)
    out.lcd3b3 <- epidish(beta.m = methDatacd3b3.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.lcd3b3, col = c("#FDD4D4", "#99E5E0"))
    out.lcd3b3
    # Choose the residual data and combined with the CpG identifiers
    cd3b4_resd <- cbind(cd3$ProbeID, as.data.frame(lm_cd3B4$residuals))
    colnames(cd3b4_resd) <- c("ProbeID", "residuals")
    rownames(cd3b4_resd) <- NULL
    cd3b4_resd_posMax <- cd3b4_resd[cd3b4_resd$residuals > 0.3,]
    cd3b4_resd_negMax <- cd3b4_resd[cd3b4_resd$residuals < -0.3,]
    cd3b4_resd_03 <- rbind(cd3b4_resd_posMax, cd3b4_resd_negMax)
    cd3b4_resd_03_data <- merge(cd3b4_resd_03, cd3, by = "ProbeID")
    cd3b4_resd_03_data_1 <- cd3b4_resd_03_data[c(1,9,10)]
    rownames(cd3b4_resd_03_data_1) <- cd3b4_resd_03_data_1[,1]
    cd3b4_resd_03_data_1 <- cd3b4_resd_03_data_1[,2:3]
    methDatacd3b4 <- cd3b4_resd_03_data_1
    methDatacd3b4.m <- as.matrix(methDatacd3b4)
    out.lcd3b4 <- epidish(beta.m = methDatacd3b4.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.lcd3b4, col = c("#FDD4D4", "#99E5E0"))
    out.lcd3b4
    # Choose the residual data and combined with the CpG identifiers
    cd3b5_resd <- cbind(cd3$ProbeID, as.data.frame(lm_cd3B5$residuals))
    colnames(cd3b5_resd) <- c("ProbeID", "residuals")
    rownames(cd3b5_resd) <- NULL
    cd3b5_resd_posMax <- cd3b5_resd[cd3b5_resd$residuals > 0.3,]
    cd3b5_resd_negMax <- cd3b5_resd[cd3b5_resd$residuals < -0.3,]
    cd3b5_resd_03 <- rbind(cd3b5_resd_posMax, cd3b5_resd_negMax)
    cd3b5_resd_03_data <- merge(cd3b5_resd_03, cd3, by = "ProbeID")
    cd3b5_resd_03_data_1 <- cd3b5_resd_03_data[c(1,11,12)]
    rownames(cd3b5_resd_03_data_1) <- cd3b5_resd_03_data_1[,1]
    cd3b5_resd_03_data_1 <- cd3b5_resd_03_data_1[,2:3]
    methDatacd3b5 <- cd3b5_resd_03_data_1
    methDatacd3b5.m <- as.matrix(methDatacd3b5)
    out.lcd3b5 <- epidish(beta.m = methDatacd3b5.m, ref.m = refData.m, method = "RPC")$estF
    boxplot(out.lcd3b5, col = c("#FDD4D4", "#99E5E0"))
    out.lcd3b5
}



