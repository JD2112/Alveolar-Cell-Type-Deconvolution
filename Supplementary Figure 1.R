# Rscript modified and re-written by Jyotirmoy Das.
# Please read the README file for the sample preparation of the analysis.

library(aplpack)
library(ggpubr)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
data <- read.table("myNorm_BAL_CD3.txt", stringsAsFactors=F, header=T, row.names = 1)
colnames(data) <- substr(colnames(data), start = 1, stop = 6)
temp <- expression(beta ~"value")
Group = factor(combineB5$grp, levels = c("SI", "BAL"))
B1_density <- ggplot(data=combineB1, aes(x=combine, group=Group, fill=Group)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum(base_family = "Times New Roman", base_size = 12) 
B11 <- B1_density + theme(legend.position = "none") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
B2_density <- ggplot(data=combineB2, aes(x=combine, group=Group, fill=Group)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum(base_family = "Times New Roman", base_size = 12)
B21 <- B2_density + theme(legend.position = "none") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
B3_density <- ggplot(data=combineB3, aes(x=combine, group=Group, fill=Group)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum(base_family = "Times New Roman", base_size = 12)
B31 <- B3_density + theme(legend.position = "none") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
B4_density <- ggplot(data=combineB4, aes(x=combine, group=Group, fill=Group)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum(base_family = "Times New Roman", base_size = 12)
B41 <- B4_density + theme(legend.position = "none") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
B5_density <- ggplot(data=combineB5, aes(x=combine, group=Group, fill=Group)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum(base_family = "Times New Roman", base_size = 12)
B51 <- B5_density + theme(legend.position = "none") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
figure1 <- ggarrange(B11, B21, B31,B41,B51,
                    labels = c("B1", "B2", "B3","B4","B5"), font.label = list(size = 12, family = "Times New Roman", face = "plain"),
                    ncol = 5, nrow = 1, common.legend = TRUE, legend = "top"
   )

# Bean plot
png("BeanPlot_mod_CD3.png", height = 6, width = 20, units = "in", res = 300)
beanplot(data,side='both', 
     border='NA', log = "", wd = 0.11,
     col = list("#FDD4D4", c("#99E5E0", "#000000")) ,
     ylab= temp,
     names = c("B1", "B2", "B3","B4","B5"),
what = c(1,1,1,0),
ll =0.04)
legend("topright", fill = c("#FDD4D4", "#99E5E0"),
    legend = c("SI", "BAL"))
dev.off()

# Bean plot
beanplot(data1,side='both', 
     border='NA', log = "", wd = 0.11,
     col = list("red", c("blue", "violet")) ,
     ylab='methylation value',
     names = c("B1", "B2", "B3","B4","B5"),	
what = c(1,1,1,0),
ll =0.04)
legend("topright", fill = c("red", "blue"),
    legend = c("IS", "BAL"))

#Density plot

library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
B1_density <- ggplot(data=combine, aes(x=combine, group=gr, fill=gr)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()
png("B1_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
B1_density
dev.off()

B2_density <- ggplot(data=combineB2, aes(x=combine, group=gr, fill=gr)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()
png("B2_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
B2_density
dev.off()

B3_density <- ggplot(data=combineB3, aes(x=combine, group=gr, fill=gr)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()
png("B3_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
B3_density
dev.off()

B4_density <- ggplot(data=combineB4, aes(x=combine, group=gr, fill=gr)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()
png("B4_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
B4_density
dev.off()

B5_density <- ggplot(data=combineB5, aes(x=combine, group=gr, fill=gr)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()
png("B5_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
B5_density
dev.off()

library(ggpubr)
figure <- ggarrange(B1_density, B2_density, B3_density,B4_density,B5_density,
                    labels = c("B1", "B2", "B3","B4","B5"),
                    ncol = 2, nrow = 3)
#figure
png("DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
figure
dev.off()
