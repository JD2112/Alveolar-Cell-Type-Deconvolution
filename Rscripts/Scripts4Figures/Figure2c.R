# Author details
## Script name: histogramPlot.R
## Purpose of the script: To visualize the beta-value distribution of two cell types mentioned in the Das et al (2021), Epigenetics
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

histogramPlot <- function(histogramPlot){
        # Load normalized data
        suppressMessages(suppressWarnings(library(ggplot2)))
        suppressMessages(suppressWarnings(library(ggpubr)))
        
        balCD3 <- read.table("myNorm_BAL_CD3.txt", 
                        stringsAsFactors = F, 
                        header = T)
        balHLADR <- read.table("myNorm_BAL_HLADR.txt", 
                        stringsAsFactors = F, 
                        header = T)

        CD3 <- as.matrix(balCD3[,2:11])# BAL data from CD3 cells
        y2 <- dnorm(CD3, mean = mean(CD3), sd = sd(CD3)) # calculate the distribution
        zz <- mean(y2) -2*sd(y2) #setting up the cut-off
        p1 <- hist(CD3)
                abline(v=zz, col = "blue", lwd =2)
                text(zz, 18 , round(zz, 1))

        HLADR <- as.matrix(balHLADR[,2:11]) # BAL data from HLA-DR cells
        y3 <- dnorm(HLADR, mean = mean(HLADR), sd = sd(HLADR)) # calculate the distribution
        zz3 <- mean(y3) -2*sd(y3) # set the cut-off
        p2 <- hist(HLADR)
                abline(v=zz3, col = "blue", lwd =2)
                text(zz3, 18 , round(zz3, 1))
        ggarrange(p1, p2, 
                common.legend = TRUE)
}

