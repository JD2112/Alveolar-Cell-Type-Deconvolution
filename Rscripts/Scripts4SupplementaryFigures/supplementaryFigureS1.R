# Author details
## Script name: supplementaryFigureS1.R
## Purpose of the script: Visualize density/bean plot from the DNA methylation deconvolution data mentioned in the Das et al (2021), Epigenetics
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

sFigure1 <- function(sFigure1){
    # Load all R bioconductor packages
    suppressMessages(suppressWarnings(library(aplpack)))
    suppressMessages(suppressWarnings(library(ggpubr)))
    suppressMessages(suppressWarnings(library(tidyverse)))
    suppressMessages(suppressWarnings(library(viridis)))
    suppressMessages(suppressWarnings(library(hrbrthemes)))
    suppressMessages(suppressWarnings(library(ggplot2)))
    suppressMessages(suppressWarnings(library(hrbrthemes)))
    suppressMessages(suppressWarnings(library(dplyr)))
    suppressMessages(suppressWarnings(library(tidyr)))
    
    data <- read.table("myNorm_BAL_CD3.txt", 
                stringsAsFactors=F, 
                header=T, 
                row.names = 1)
    colnames(data) <- substr(colnames(data), start = 1, stop = 6)
    temp <- expression(beta ~"value")
    Group = factor(combineB5$grp, levels = c("SI", "BAL"))
    # Density plot for Smaple P1
    B1_density <- ggplot(data=combineB1, aes(x=combine, group=Group, fill=Group)) +
                    geom_density(adjust=1.5, alpha=.4) +
                    theme_ipsum(base_family = "Times New Roman", 
                        base_size = 12) 
    B11 <- B1_density + 
        theme(legend.position = "none") + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    # Density plot for Smaple P2
    B2_density <- ggplot(data=combineB2, aes(x=combine, group=Group, fill=Group)) +
                    geom_density(adjust=1.5, alpha=.4) +
                    theme_ipsum(base_family = "Times New Roman", 
                        base_size = 12)
    B21 <- B2_density + theme(legend.position = "none") + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    # Density plot for Smaple P3
    B3_density <- ggplot(data=combineB3, aes(x=combine, group=Group, fill=Group)) +
                    geom_density(adjust=1.5, alpha=.4) +
                    theme_ipsum(base_family = "Times New Roman", base_size = 12)
    B31 <- B3_density + theme(legend.position = "none") + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    # Density plot for Smaple P4
    B4_density <- ggplot(data=combineB4, aes(x=combine, group=Group, fill=Group)) +
                    geom_density(adjust=1.5, alpha=.4) +
                    theme_ipsum(base_family = "Times New Roman", 
                        base_size = 12)
    B41 <- B4_density + theme(legend.position = "none") + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    # Density plot for Smaple P5
    B5_density <- ggplot(data=combineB5, aes(x=combine, group=Group, fill=Group)) +
                    geom_density(adjust=1.5, alpha=.4) +
                    theme_ipsum(base_family = "Times New Roman", 
                        base_size = 12)
    B51 <- B5_density + theme(legend.position = "none") + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    # arrange all figures together
    figure1 <- ggarrange(B11, B21, B31,B41,B51,
                    labels = c("B1", "B2", "B3","B4","B5"), 
                    font.label = list(size = 12, 
                        family = "Times New Roman", 
                            face = "plain"),
                    ncol = 5, 
                    nrow = 1, 
                    common.legend = TRUE, legend = "top"
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
    #Density plot for Sample P1
    B1_density <- ggplot(data=combine, aes(x=combine, group=gr, fill=gr)) +
            geom_density(adjust=1.5, alpha=.4) +
            theme_ipsum()
    png("B1_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
    B1_density
    dev.off()
    #Density plot for Sample P2
    B2_density <- ggplot(data=combineB2, aes(x=combine, group=gr, fill=gr)) +
            geom_density(adjust=1.5, alpha=.4) +
            theme_ipsum()
    png("B2_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
    B2_density
    dev.off()
    #Density plot for Sample P3
    B3_density <- ggplot(data=combineB3, aes(x=combine, group=gr, fill=gr)) +
            geom_density(adjust=1.5, alpha=.4) +
            theme_ipsum()
    png("B3_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
    B3_density
    dev.off()
    #Density plot for Sample P4
    B4_density <- ggplot(data=combineB4, aes(x=combine, group=gr, fill=gr)) +
            geom_density(adjust=1.5, alpha=.4) +
            theme_ipsum()
    png("B4_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
    B4_density
    dev.off()
    #Density plot for Sample P5
    B5_density <- ggplot(data=combineB5, aes(x=combine, group=gr, fill=gr)) +
            geom_density(adjust=1.5, alpha=.4) +
            theme_ipsum()
    png("B5_DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
    B5_density
    dev.off()
    figure <- ggarrange(B1_density, B2_density, B3_density,B4_density,B5_density,
                    labels = c("B1", "B2", "B3","B4","B5"),
                    ncol = 2, nrow = 3)
    #save the figure
    png("DensityPlot_CD3.png", height = 15, width = 10, units = "in", res = 300)
    figure
    dev.off()
}