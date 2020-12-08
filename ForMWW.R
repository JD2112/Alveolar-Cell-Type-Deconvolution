# R scripts for hypothesis testing
# New analysis: Mann-Whitney-Wilcoxon test for hypothsis
# H0: there is no difference in the BAL vs SI samples over all participants
# Ha: there is a difference in the BAL vs SI samples over all participants

# Procedure:
# Combine all participants all CpGs beta values separately for HLA-DR and CD3 cell types
setwd("/home/jyoda68/Desktop/Link to BackUp_All/BALdata/BALdata/Part2/DataSheets/")
hladrCPGs <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors = F, header = T)
head(hladrCPGs)
# Separate si and bal data
hladr_si <- hladrCPGs[, c(1,2,4,6,8,10)]
hladr_bal <- hladrCPGs[, c(1,3,5,7,9,11)]
# avarge si and bal data
hladr_si$avg <- (hladr_si$B1_SPU_HLDADR+hladr_si$B2_SPU_HLDADR + hladr_si$B3_SPU_HLDADR + 
                   hladr_si$B4_SPU_HLDADR + hladr_si$B5_SPU_HLDADR)/5
head(hladr_si)
hladr_bal$avg <- (hladr_bal$B1_LAV_HLDADR+hladr_bal$B2_LAV_HLDADR + hladr_bal$B3_LAV_HLDADR + 
                   hladr_bal$B4_LAV_HLDADR + hladr_bal$B5_LAV_HLDADR)/5
head(hladr_bal)
# Create a second file for MWW test
hladr_siAVG <- hladr_si[,c(1,7)]
hladr_siAVG$type <- rep("SI")
head(hladr_siAVG)
hladr_balAVG <- hladr_bal[,c(1,7)]
hladr_balAVG$type <- rep("BAL")
head(hladr_balAVG)

# combine two final files 
hladr_si_bal <- rbind(hladr_siAVG, hladr_balAVG)
head(hladr_si_bal)
tail(hladr_si_bal)
# Perform MWW test
wilcox.test(avg ~ type, hladr_si_bal)
kruskal.test(avg ~ type, hladr_si_bal)

# plotting the data
library(ggplot2)
library(RColorBrewer)
# Basic violin plot
p <- ggplot(hladr_si_bal, aes(x=type, y=avg)) + 
  geom_violin()
p + geom_boxplot(width=0.1)
p + stat_summary(fun.data=mean_sdl, mult=1, 
                 geom="pointrange", color="red")

# Violin plot (boxplot inset)
p <- ggplot(hladr_si_bal, aes(x=type, y=avg, fill = type)) + 
  geom_violin(trim=FALSE) + #fill=c("#614051", "2A52BE"))+
  labs(x="", y = "beta value")+
  geom_boxplot(width=0.01)+
  #geom_hline(yintercept = mean(hladr_si_bal$avg), linetype = 2)+
  theme_classic()
p + scale_fill_manual(values = c("#FDD4D4", "#99E5E0")) +
  theme(axis.text.x = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.text.y = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(beta~"value")) +
  ggtitle("a") + 
  theme(text = element_text(size = 12, family = "Times New Roman", color = "black")) +
  stat_compare_means(label = "p.signif", paired = T, label.x = 1.5)+
  geom_segment(aes(x= 1.2, y = 0.98, xend = 1.8, yend = 0.98))
# Statistical calculation for adding the significance level
mycom <- compare_means(avg ~ type,  data = hladr_si_bal)
install.packages("ggpubr")
library(ggpubr)

# Violin plot with dots
p <- ggplot(hladr_si_bal, aes(x=type, y=avg, fill = type)) + 
  geom_violin(trim=FALSE) + #fill=c("#614051", "2A52BE"))+
  labs(x="", y = "beta value")+
  #geom_boxplot(width=0.01)+
  #geom_hline(yintercept = mean(hladr_si_bal$avg), linetype = 2)+
  theme_classic()
p2 <- p + scale_fill_manual(values = c("#FDD4D4", "#99E5E0")) +
  theme(axis.text.x = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.text.y = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(beta~"value")) +
  ggtitle("a (HLA-DR)") + 
  theme(text = element_text(size = 12, family = "Times New Roman", color = "black")) +
  stat_compare_means(label = "p.signif", paired = T, label.x = 1.5)+
  geom_segment(aes(x= 1.2, y = 0.98, xend = 1.8, yend = 0.98)) +
  #stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1)) +
  stat_summary(fun.data=data_summary, show.legend = T) +
  stat_summary(fun.y = mean, color = "black", geom ="point", aes(group = 1), size = 5,
               show.legend = FALSE)

foo <- as.data.frame(ggplot_build(p2)$data[[4]])
p2 +
  annotate("label", x = foo$x + 0.14, y = 0.47, color = "black", label = round(foo$y, digits = 4))
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


# CD3 data
CD3CPGs <- read.table("myNorm_BAL_CD3.txt", stringsAsFactors = F, header = T)
head(CD3CPGs)
# Separate si and bal data
cd3_si <- CD3CPGs[, c(1,2,4,6,8,10)]
cd3_bal <- CD3CPGs[, c(1,3,5,7,9,11)]
# avarge si and bal data
cd3_si$avg <- (cd3_si$B1_SPU_CD3+cd3_si$B2_SPU_CD3 + cd3_si$B3_SPU_CD3 + 
                   cd3_si$B4_SPU_CD3 + cd3_si$B5_SPU_CD3)/5
head(cd3_si)
cd3_bal$avg <- (cd3_bal$B1_LAV_CD3+cd3_bal$B2_LAV_CD3 + cd3_bal$B3_LAV_CD3 + 
                    cd3_bal$B4_LAV_CD3 + cd3_bal$B5_LAV_CD3)/5
head(cd3_bal)
# Create a second file for MWW test
cd3_siAVG <- cd3_si[,c(1,7)]
cd3_siAVG$type <- rep("SI")
head(cd3_siAVG)
cd3_balAVG <- cd3_bal[,c(1,7)]
cd3_balAVG$type <- rep("BAL")
head(cd3_balAVG)

# combine two final files 
cd3_si_bal <- rbind(cd3_siAVG, cd3_balAVG)
head(cd3_si_bal)
tail(cd3_si_bal)
# Perform MWW test
wilcox.test(avg ~ type, cd3_si_bal)
kruskal.test(avg ~ type, cd3_si_bal)

# Violin plot (boxplot inset)
p <- ggplot(cd3_si_bal, aes(x=type, y=avg, fill = type)) + 
  geom_violin(trim=FALSE) + #fill=c("#614051", "2A52BE"))+
  labs(x="", y = "beta value")+
  geom_boxplot(width=0.01)+
  #geom_hline(yintercept = mean(hladr_si_bal$avg), linetype = 2)+
  theme_classic()
p + scale_fill_manual(values = c("#FDD4D4", "#99E5E0")) +
  theme(axis.text.x = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.text.y = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(beta~"value")) +
  ggtitle("b") + 
  theme(text = element_text(size = 12, family = "Times New Roman", color = "black")) +
  stat_compare_means(label = "p.signif", paired = T, label.x = 1.5)+
  geom_segment(aes(x= 1.2, y = 0.98, xend = 1.8, yend = 0.98))
# Statistical calculation for adding the significance level
mycom <- compare_means(avg ~ type,  data = cd3_si_bal)
#install.packages("ggpubr")
library(ggpubr)

# Violin plot with dots
p <- ggplot(cd3_si_bal, aes(x=type, y=avg, fill = type)) + 
  geom_violin(trim=FALSE) + #fill=c("#614051", "2A52BE"))+
  labs(x="", y = "beta value")+
  #geom_boxplot(width=0.01)+
  #geom_hline(yintercept = mean(hladr_si_bal$avg), linetype = 2)+
  theme_classic()
p2 <- p + scale_fill_manual(values = c("#FDD4D4", "#99E5E0")) +
  theme(axis.text.x = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.text.y = element_text(family = "Times New Roman", size = 12, colour = "black"))+
  theme(axis.title.y = element_text(family = "Times New Roman", size = 14, colour = "black", face = "bold"))+
  theme(legend.position = "none")+
  scale_y_continuous(name = expression(beta~"value")) +
  ggtitle("b (CD3)") + 
  theme(text = element_text(size = 12, family = "Times New Roman", color = "black")) +
  stat_compare_means(label = "p.signif", paired = T, label.x = 1.5)+
  geom_segment(aes(x= 1.2, y = 0.98, xend = 1.8, yend = 0.98)) +
  #stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1)) +
  stat_summary(fun.data=data_summary, show.legend = T) +
  stat_summary(fun.y = mean, color = "black", geom ="point", aes(group = 1), size = 5,
               show.legend = FALSE)

foo <- as.data.frame(ggplot_build(p2)$data[[4]])
p2 +
  annotate("label", x = foo$x + 0.14, y = 0.47, color = "black", label = round(foo$y, digits = 4))
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Re-analysis of Ringh et al data with average beta values
# Training data = average of all participants 2763 CpG beta values
# Test data = 850K data set from Ringh et al
library(EpiDISH)
# data(centEpiFibIC.m)
# data(DummyBeta.m)
referenceData <- read.table("NewTrainData.txt", stringsAsFactors = F, header = T, row.names = 1)
methData <- read.table("Stockholm_normalizedValue_GSE133062.txt", stringsAsFactors = F, header = T, row.names = 1)
methData.m <- as.matrix(methData)
refData.m <- as.matrix(referenceData)
out.lB1 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
# Plot the graph
install.packages('extrafont')
library(extrafont)
pdf("Figure 4 validation.pdf", family="Times New Roman", width=3, height=4)
outer = FALSE
line = 1
cex = 1.3
adj  = 0.025
par(mfrow=c(1,1), tcl=-0.5, family="Times New Roman", omi=c(0.3,0.3,0,0), 
    mar=c(2.3, 2.3, 1.8, 1.8), oma = c(2, 2, 1, 1))
boxplot(out.lB1, col = c("#FDD4D4", "#99E5E0"), font.family = "Times New Roman")
mtext("Cell proportion", side = 2, outer = T, at = 0.5, cex = 1.3, family = "Times New Roman")
mtext("Cell types", side=1, outer=T, at=0.5, cex = 1.3, family = "Times New Roman")
dev.off()

### Analysis with Cystic Fibrosis data
setwd("../../Part3-ValidationWithCysticFibrosis/")
CFdata <- read.table("normalizedValue_GSE132547.txt", stringsAsFactors = F, header = T)
  