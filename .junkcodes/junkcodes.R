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


# Violin Plot
library(tidyverse)
library(viridis)
combine %>%
  ggplot( aes(x=gr, y=combine, fill=gr)) +
    geom_violin() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("Violin chart") +
    xlab("")

data1 <- data %>% 
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = round(as.numeric(value),0)) %>%
  filter(text %in% c("B1_SPU","B1_LAV","B2_SPU","B2_LAV","B3_SPU", "B3_LAV", "B4_SPU", "B4_LAV", "B5_SPU","B5_LAV"))

# Plot
p <- data1 %>%
  mutate(text = fct_reorder(text, value)) %>% # Reorder data
  ggplot( aes(x=text, y=value, fill=text, color=text)) +
    geom_violin(width=2.1, size=0.2) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none"
    ) +
    coord_flip() + # This switch X and Y axis and allows to get the horizontal version
    xlab("") +
    ylab("")

p



# Calculating the wilcox test (two.sided)

comH <- data.frame(c(boxH[,"SI"], boxH[, "BAL"]))
colnames(comH) <- "combine"
#combineB1$gr <- factor(rep(1:2, each = nrow(data)))
comH$grp <- factor(rep(c("SI", "BAL"), each = nrow(boxH)))
wilcox.test(combine ~ grp, data = comH)

comC <- data.frame(c(boxC[,"SI"], boxC[, "BAL"]))
colnames(comC) <- "combine"
comC$grp <- factor(rep(c("SI", "BAL"), each = nrow(boxC)))
wilcox.test(combine ~ grp, data = comC)

> t.test(combine ~ grp, data = comH)


> t.test(combine ~ grp, data = comC)


# Heatmap
library(ComplexHeatmap)
hm <- as.matrix(combine)

colnames(hm) <- c("", "HLA-DR", "", "CD3")
Heatmap(hm, column_title_gp = gpar(fontsize = 12, fontfamily = "Times New Roman"),column_names_rot = 0, row_title_gp = gpar(fontsize = 12, fontfamily = "Times New Roman"), heatmap_legend_param = list(title = "est", color_bar = "discrete"))




### Analysis with Cystic Fibrosis data
setwd("../../Part3-ValidationWithCysticFibrosis/")
CFdata <- read.table("normalizedValue_GSE132547.txt", stringsAsFactors = F, header = T)

  #install.packages("PerformanceAnalytics")
  # load packages
  suppressMessages(suppressWarnings(library("PerformanceAnalytics")))
  suppressMessages(suppressWarnings(library("corrplot")))
  suppressMessages(suppressWarnings(library("ggpubr")))
  suppressMessages(suppressWarnings(library("ggplot2")))
  # correlation chart
  chart.Correlation(data1, histogram=TRUE, pch=19)
  # correlation plot
  corrplot(data, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
  # correaltion coefficient calculation
  cor(x, y, method = c("pearson", "kendall", "spearman"))
  cor.test(x, y, method=c("pearson", "kendall", "spearman"))
  # scatter plot 
  # scatter plot for P1
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
corrplot.mixed(res2$r, 
            color="gray", 
            fill = "gray", 
            order="hclust",
            number.cex=1, 
            tl.cex = 1, 
            p.mat = res2$P, 
            sig.level = 1e-16, 
            insig = "blank")

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

data2  <- cor(data, method = "s")
corrplot(data2, method = "number", order = "hclust", 
         tl.col = "black", tl.srt = 45)
