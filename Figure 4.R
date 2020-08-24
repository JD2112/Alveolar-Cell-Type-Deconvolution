# Rscript modified and re-written by Jyotirmoy Das.
# Please read the README file for the sample preparation of the analysis.

setwd("~/Desktop/Link to BackUp_All/BALdata/BALdata/Part2/DataSheets")
uniqCpGs <- read.table("cell-specific-CpGs_ourStudy.txt", sstringsAsFactors = F, header = T)
uniqCpGs <- read.table("cell-specific-CpGs_ourStudy.txt", stringsAsFactors = F, header = T)
head(uniqCpGs)
BAL_HLADR <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors = F, header = T)
BAL_CD3 <- read.table("myNorm_BAL_CD3.txt", stringsAsFactors = F, header = T)
uniqCpGsHLADR <- uniqCpGs[1:594,]
uniqCpGsCD3 <- uniqCpGs[595:2886,]

library(EpiDISH)
# data(centEpiFibIC.m)
# data(DummyBeta.m)
referenceData <- read.table("refDataB1_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
methData <- read.table("Stockholm_normalizedValue_GSE133062.txt", stringsAsFactors = F, header = T, row.names = 1)
out.l <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
methData.m <- as.matrix(methData)
refData.m <- as.matrix(referenceData)
out.lB1 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
boxplot(out.l)
boxplot(out.lB1, col = c("#FDD4D4", "#99E5E0"))

referenceData <- read.table("refDataB2_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
refData.m <- as.matrix(referenceData)
out.lB2 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
boxplot(out.lB2, col = c("#FDD4D4", "#99E5E0"))

referenceData <- read.table("refDataB3_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
refData.m <- as.matrix(referenceData)
out.lB3 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
boxplot(out.lB3, col = c("#FDD4D4", "#99E5E0"))

referenceData <- read.table("refDataB4_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
refData.m <- as.matrix(referenceData)
out.lB4 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
boxplot(out.lB4, col = c("#FDD4D4", "#99E5E0"))

referenceData <- read.table("refDataB5_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
refData.m <- as.matrix(referenceData)
out.lB5 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
boxplot(out.lB5, col = c("#FDD4D4", "#99E5E0"))

referenceData <- read.table("refDataAVG_CD3_HLADR.txt", stringsAsFactors = F, header = T, row.names = 1)
refData.m <- as.matrix(referenceData)
out.lBAVG <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
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
