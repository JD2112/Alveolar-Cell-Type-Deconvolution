library("Hmisc")
library("ggpubr")
data <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors=F, header=T)
colnames(data) <- substr(colnames(data), start = 1, stop = 6)
dataCD3 <- read.table("~/../CD3/myNorm_BAL_CD3.txt", stringsAsFactors=F, header=T)
colnames(dataCD3) <- substr(colnames(dataCD3), start = 1, stop = 6)
# Using ggscatter only
B1 <- ggscatter(data, x = "B1_SPU", y = "B1_LAV", xlab = FALSE, ylab = FALSE,
                color = "gray", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "black"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                conf.int.level = 0.95,
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                cor.coef.size = 6
) 

B2 <- ggscatter(data, x = "B2_SPU", y = "B2_LAV", xlab = FALSE, ylab = FALSE,
                color = "gray", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                conf.int.level = 0.95,
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                cor.coef.size = 6
)

B3 <- ggscatter(data, x = "B3_SPU", y = "B3_LAV", xlab = FALSE, ylab = FALSE,
                color = "gray", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                conf.int.level = 0.95,
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                cor.coef.size = 6
) 
B4 <- ggscatter(data, x = "B4_SPU", y = "B4_LAV", xlab = FALSE, ylab = FALSE,
                color = "gray", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                conf.int.level = 0.95,
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                cor.coef.size = 6
) 

B5 <- ggscatter(data, x = "B5_SPU", y = "B5_LAV", xlab = FALSE, ylab = FALSE,
                color = "gray", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                conf.int.level = 0.95,
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                cor.coef.size = 6
) 

B1_CD3 <- ggscatter(dataCD3, x = "B1_SPU", y = "B1_LAV", xlab = FALSE, ylab = FALSE,
                    color = "gray", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    conf.int.level = 0.95,
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                    cor.coef.size = 6
) 
B2_CD3 <- ggscatter(dataCD3, x = "B2_SPU", y = "B2_LAV", xlab = FALSE, ylab = FALSE,
                    color = "gray", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    conf.int.level = 0.95,
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                    cor.coef.size = 6
) 
B3_CD3 <- ggscatter(dataCD3, x = "B3_SPU", y = "B3_LAV", xlab = FALSE, ylab = FALSE,
                    color = "gray", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    conf.int.level = 0.95,
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                    cor.coef.size = 6
) 

B4_CD3 <- ggscatter(dataCD3, x = "B4_SPU", y = "B4_LAV", xlab = FALSE, ylab = FALSE,
                    color = "gray", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    conf.int.level = 0.95,
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                    cor.coef.size = 6
) 


B5_CD3 <- ggscatter(dataCD3, x = "B5_SPU", y = "B5_LAV", xlab = FALSE, ylab = FALSE,
                    color = "gray", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    conf.int.level = 0.95,
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "spearman", label.sep = "\n", color = "black"),
                    cor.coef.size = 6
) 

figure2 <- ggarrange(B1, B2, B3, B4, B5, B1_CD3, B2_CD3, B3_CD3, B4_CD3, B5_CD3, 
                     labels = c("B1", "B2", "B3", "B4", "B5", "B1", "B2", "B3", "B4", "B5"),
                     hjust = -2.5,
                     vjust = 2,
                     ncol = 5, nrow = 2)
  annotate_figure(figure2,
                  left = text_grob("Bronchoalveolar Lavage (BAL)", color = "black", family = "Times New Roman", size = 16, rot = 90),
                  bottom = text_grob("Sputum Induction (SI)", color = "black", family = "Times New Roman", size = 16))
  
  
  # Using stat_cor() with ggscatter
  B1 <- ggscatter(data, x = "B1_SPU", y = "B1_LAV",xlab = FALSE, ylab = FALSE,
                color = "black", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "black"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                #cor.coeff.args = list(method = "spearman")
)
B1_1 <- B1 + stat_cor(aes(color = "#0000FF", 
                          label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                          family = "Times New Roman"), 
                      size = 6,
                      show.legend = F)
B2 <- ggscatter(data, x = "B2_SPU", y = "B2_LAV",xlab = FALSE, ylab = FALSE,
                color = "black", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                #cor.coeff.args = list(method = "spearman")
)
B2_1 <- B2 + stat_cor(aes(color = "#0000FF", 
                          label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                          family = "Times New Roman"), 
                      size = 6,
                      show.legend = F)
B3 <- ggscatter(data, x = "B3_SPU", y = "B3_LAV",xlab = FALSE, ylab = FALSE,
                color = "black", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                #cor.coeff.args = list(method = "spearman")
)
B3_1 <- B3 + stat_cor(aes(color = "#0000FF", 
                          label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                          family = "Times New Roman"), 
                      size = 6,
                      show.legend = F)
B4 <- ggscatter(data, x = "B4_SPU", y = "B4_LAV",xlab = FALSE, ylab = FALSE,
                color = "black", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                #cor.coeff.args = list(method = "spearman")
)
B4_1 <- B4 + stat_cor(aes(color = "#0000FF", 
                          label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                          family = "Times New Roman"), 
                      size = 6,
                      show.legend = F)

B5 <- ggscatter(data, x = "B5_SPU", y = "B5_LAV",xlab = FALSE, ylab = FALSE,
                color = "black", shape = 19,  # Points color, shape and size
                add = "reg.line",  # Add regressin line,
                font.family = "Times New Roman", font.label =12,
                add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                #cor.coeff.args = list(method = "spearman")
)
B5_1 <- B5 + stat_cor(aes(color = "#0000FF", 
                          label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                          family = "Times New Roman"), 
                      size = 6,
                      show.legend = F) 
B1_CD3 <- ggscatter(dataCD3, x = "B1_SPU", y = "B1_LAV",xlab = FALSE, ylab = FALSE,
                    color = "black", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #cor.coeff.args = list(method = "spearman")
)
B1_CD3_1 <- B1_CD3 + stat_cor(aes(color = "#0000FF", 
                                  label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                                  family = "Times New Roman"), 
                              size = 6,
                              show.legend = F) 

B2_CD3 <- ggscatter(dataCD3, x = "B2_SPU", y = "B2_LAV",xlab = FALSE, ylab = FALSE,
                    color = "black", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #cor.coeff.args = list(method = "spearman")
)
B2_CD3_1 <- B2_CD3 + stat_cor(aes(color = "#0000FF", 
                                  label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                                  family = "Times New Roman"), 
                              size = 6,
                              show.legend = F) 

B3_CD3 <- ggscatter(dataCD3, x = "B3_SPU", y = "B3_LAV",xlab = FALSE, ylab = FALSE,
                    color = "black", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #cor.coeff.args = list(method = "spearman")
)

B3_CD3_1 <- B3_CD3 + stat_cor(aes(color = "#0000FF", 
                                  label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                                  family = "Times New Roman"), 
                              size = 6,
                              show.legend = F) 

B4_CD3 <- ggscatter(dataCD3, x = "B4_SPU", y = "B4_LAV",xlab = FALSE, ylab = FALSE,
                    color = "black", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "lightblue"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #cor.coeff.args = list(method = "spearman")
)

B4_CD3_1 <- B4_CD3 + stat_cor(aes(color = "#0000FF", 
                                  label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                                  family = "Times New Roman"), 
                              size = 6,
                              show.legend = F) 
B5_CD3 <- ggscatter(dataCD3, x = "B5_SPU", y = "B5_LAV", xlab = FALSE, ylab = FALSE,
                    color = "black", shape = 19,  # Points color, shape and size
                    add = "reg.line",  # Add regressin line,
                    font.family = "Times New Roman", font.label =12,
                    add.params = list(color = "blue", fill = "gray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    #conf.int.level = 0.95,
                    #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    #cor.coeff.args = list(method = "spearman"), #label.sep = "\n", color = "blue"),
                    #cor.coef.size = 14
) 
B5_CD3_1 <- B5_CD3 + stat_cor(aes(color = "#0000FF", 
                                  label = paste(..rr.label.., ..p.label..,sep = "~`;`~"), 
                                  family = "Times New Roman"), 
                              size = 6,
                              show.legend = F) 

figure2 <- ggarrange(B1_1, B2_1, B3_1, B4_1, B5_1, B1_CD3_1, B2_CD3_1, B3_CD3_1, B4_CD3_1, B5_CD3_1 + rremove("x.text") + rremove("y.text"), 
                     labels = c("B1", "B2", "B3", "B4", "B5", "B1", "B2", "B3", "B4", "B5"),
                     ncol = 5, nrow = 2)
annotate_figure(figure2,
                left = text_grob("Bronchoalveolar Lavage (BAL)", color = "black", family = "Times New Roman", size = 14, rot = 90),
                bottom = text_grob("Sputum Induction (SI)", color = "black", family = "Times New Roman", size = 14))


## Alternative to regression plot using the caret package with linear regression model analysis
library(tidyverse)
library(caret)

data <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors=F, header = T)

head(data)

p1 <- ggplot(data, aes(x = B1_SPU_HLDADR, y = B1_LAV_HLDADR)) +
  geom_point() +
  stat_smooth()
p1 + theme(text = element_text(family = "Times New Roman", face = "bold", size = 14), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), axis.line = element_line(size = 1),
        panel.grid.major =  element_blank(), panel.grid.minor = element_blank(),
  ) +
  labs(x = "IS", y = "BAL")


model <- lm(B1_SPU_HLDADR ~ B1_LAV_HLDADR, data = data)
summary(model)$coef

summary(model)

png("B5_HLA-DR_BALvsIS.png", height = 10, width = 8, unit = "in", res = 300)
p1 <- ggplot(data, aes(x = B5_SPU_HLDADR, y = B5_LAV_HLDADR)) +
  geom_point() +
  stat_smooth()
p1 + theme(text = element_text(family = "Times New Roman", face = "bold", size = 14), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), axis.line = element_line(size = 1),
        panel.grid.major =  element_blank(), panel.grid.minor = element_blank(),
  ) +
  labs(x = "IS", y = "BAL")
dev.off()

o6 <- ggplot(data, aes(x = B5_SPU_HLDADR, y = B5_LAV_HLDADR)) +
  geom_point(alpha = 0.1) +
  geom_rug(alpha = 0.01) +
	stat_smooth()

library(ggplot2)
# install.packages("ggpointdensity")
library(ggpointdensity)
png("B5_HLA-DR_BALvsIS-viridis1.png", height = 10, width = 8, unit = "in", res = 300)
p2 <- ggplot(data, aes(x = B5_SPU_HLDADR, y = B5_LAV_HLDADR)) + 
	geom_pointdensity() + scale_color_viridis_c() +
	stat_smooth()
p2 + theme(text = element_text(family = "Times New Roman", face = "bold", size = 14), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), axis.line = element_line(size = 1),
        panel.grid.major =  element_blank(), panel.grid.minor = element_blank(),
  ) +
  labs(x = "IS", y = "BAL")
dev.off()

library(ggplot2)
library(ggsubplot)

# Scatterplot with subplots (simple)
ggplot(dat, aes(x=xvar, y=yvar)) +
  geom_point(shape=1) +
  geom_subplot2d(aes(xvar, yvar,
                     subplot = geom_bar(aes(rep("dummy", length(xvar)), ..count..))), bins = c(15,15), ref = NULL, width = rel(0.8), ply.aes = FALSE)



# Final Graph calculation
ggplotRegression <- function (fit) {

require(ggplot2)

lb1 <- paste("R^2 == ", signif(summary(fit)$adj.r.squared, 4))
lb2 <- paste("Intercept == ", signif(fit$coef[[1]],4))
lb3 <- paste("Slope == ", signif(fit$coef[[2]], 4))
lb4 <- paste(expression(italic('p')),"-value == ", 6.35e-63)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 14), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12), axis.line = element_line(size = 1),
        panel.grid.major =  element_blank(), panel.grid.minor = element_blank(),
  ) +
  labs(x = "IS", y = "BAL") +
annotate("label", x= 0.7, y = 0.19, label = lb1, parse =TRUE, hjust =0,vjust = 0.5,label.size = NA, label.padding = unit (0.5, "lines"),label.r = unit(0.05, "lines")) +
annotate("label", x= 0.7, y = 0.15, label = lb2, parse =TRUE, hjust =0, label.size = NA,vjust = 0.5 , label.padding = unit (0.5, "lines"),label.r = unit(0.05, "lines")) +
annotate("label", x= 0.7, y = 0.11, label = lb3,parse =TRUE, hjust =0, label.size = NA,vjust = 0.5 , label.padding = unit (0.5, "lines"), label.r = unit(0.05, "lines"))+
annotate("label", x= 0.7, y = 0.07, label = lb4, fontface = "italic", parse =TRUE, hjust =0, vjust = 0.5 , label.padding = unit (0.5, "lines"),label.size = NA, label.r = unit(0.05, "lines"))
}

fit1 <- lm(B1_SPU_CD3 ~ B1_LAV_CD3, data = data)
png("B1_CD3_BALvsIS-1.png", height = 10, width = 8, unit = "in", res = 300)
ggplotRegression(fit1)
dev.off()

fit2 <- lm(B2_SPU_CD3 ~ B2_LAV_CD3, data = data)
png("B2_CD3_BALvsIS-1.png", height = 10, width = 8, unit = "in", res = 300)
ggplotRegression(fit2)
dev.off()

fit3 <- lm(B3_SPU_CD3 ~ B3_LAV_CD3, data = data)
png("B3_CD3_BALvsIS-1.png", height = 10, width = 8, unit = "in", res = 300)
ggplotRegression(fit3)
dev.off()

fit4 <- lm(B4_SPU_CD3 ~ B4_LAV_CD3, data = data)
png("B4_CD3_BALvsIS-1.png", height = 10, width = 8, unit = "in", res = 300)
ggplotRegression(fit4)
dev.off()

fit5 <- lm(B5_SPU_CD3 ~ B5_LAV_CD3, data = data)
png("B5_CD3_BALvsIS-1.png", height = 10, width = 8, unit = "in", res = 300)
ggplotRegression(fit5)
dev.off()
