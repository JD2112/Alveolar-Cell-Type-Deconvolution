# Author details
## Script name: supplementaryFigureS2.R
## Purpose of the script: Correlation plot from the DNA methylation deconvolution data mentioned in the Das et al (2021), Epigenetics
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

# Rscript modified and re-written by Jyotirmoy Das.
# Please read the README file for the sample preparation of the analysis.

sFigure2 <- function(sFigure2){
  suppressMessages(suppressWarnings(library("Hmisc")))
  suppressMessages(suppressWarnings(library(corrplot)))
  res2 <- rcorr(as.matrix(data1), type = "pearson") # type = c("spearman", "pearson")
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
}