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
