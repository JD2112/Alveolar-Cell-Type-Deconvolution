# Rscript modified and re-written by Jyotirmoy Das.
# Please read the README file for the sample preparation of the analysis.

install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")

chart.Correlation(data1, histogram=TRUE, pch=19)

library(corrplot)
corrplot(data, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

library("ggpubr")
cor(x, y, method = c("pearson", "kendall", "spearman"))
cor.test(x, y, method=c("pearson", "kendall", "spearman"))

ggscatter(data, x = "B1_SPU", y = "B1_LAV", color="gray", fill = "gray",
          add = "reg.line", conf.int = TRUE, 
	  add.params = list(color = "blue", fill = "black"), 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Induced Sputum", ylab = "BAL")

B1 <- ggscatter(data, x = "B1_SPU", y = "B1_LAV",
   color = "gray", shape = 19,  # Points color, shape and size
   add = "reg.line",  # Add regressin line,
	font.family = "Times New Roman", font.label =12,
   add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman")
   )

png("B1_CD3_BALvsIS-2.png", height = 10, width = 8, unit = "in", res = 300)
B1
dev.off()

ad.test(data$B1_SPU)
ad.test(data$B1_LAV)
ad.test(data$B2_SPU)
ad.test(data$B2_LAV)
ad.test(data$B3_SPU)
ad.test(data$B3_LAV)
ad.test(data$B4_SPU)
ad.test(data$B4_LAV)
ad.test(data$B5_SPU)
ad.test(data$B5_LAV)

ggscatter(data, x = "B1_SPU", y = "B1_LAV",
   add = "loess", conf.int = TRUE)

B5 <- cor.test(data$B5_SPU, data$B5_LAV, 
                    method = "spearman")
B5

library("Hmisc")
res2 <- rcorr(as.matrix(data))
res2

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

flattenCorrMatrix(res2$r, res2$P)

library(corrplot)
corrplot.mixed(res2$r, color="gray", fill = "gray", order="hclust",number.cex=1, tl.cex = 1, p.mat = res2$P, sig.level = 1e-16, insig = "blank")

# Updated 2021-02-21
load("~/Documents/ML/BALvsIS/BALdata_Last/Part2_ValidationWithStockholmData/DataSheets/.RData")
setwd("~/Documents/ML/BALvsIS/BALdata_Last/Part2_ValidationWithStockholmData/DataSheets")
head(data)
library(corrplot)

data1 <- data
colnames(data1) <- gsub(x = colnames(data1), pattern = "B1_", replacement = "P1_")
colnames(data1) <- gsub(x = colnames(data1), pattern = "B2_", replacement = "P2_")
colnames(data1) <- gsub(x = colnames(data1), pattern = "B3_", replacement = "P3_")
colnames(data1) <- gsub(x = colnames(data1), pattern = "B4_", replacement = "P4_")
colnames(data1) <- gsub(x = colnames(data1), pattern = "B5_", replacement = "P5_")

data2  <- cor(data, method = "s)
corrplot(data2, method = "number", order = "hclust", 
         tl.col = "black", tl.srt = 45)


library("Hmisc")
res2 <- rcorr(as.matrix(data1))
res2

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flattenCorrMatrix(res2$r, res2$P)


library(corrplot)
par(family = "Times New Roman")
corrplot(res2$r, color="gray", 
               fill = "gray",
               method = "number",
               order="hclust",
               number.cex=1, 
               tl.cex = 1, 
               p.mat = res2$P, 
              tl.col = "black",
               sig.level = 1e-16, insig = "blank",
              tl.srt = 45)

