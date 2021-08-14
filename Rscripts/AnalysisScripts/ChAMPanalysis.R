# Author details
## Script name: ChAMPanalysis.R
## Purpose of the script: To analyze the DNA methylation data mentioned in the Das et al (2021), Epigenetics
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

ChAMPanalysis450K <- function(ChAMPanalysis450K){
    suppressMessages(suppressWarnings(library(ChAMP)))
    myLoad <- champ.load()
    myNorm <- champ.norm(beta=myLoad$beta, 
            arraytype = "450K", 
            cores = 30) # Use BMIQ normalization
    myDMP <- champ.DMP(beta = myNorm, 
            pheno=myLoad$pd$Sample_Group) # use BH-corrected p-value < 0.05
    # Save normalized beta values
    write.table(myNorm, file = "myNorm_BAL_HLADR.txt", 
            sep = "\t", 
            quote = FALSE)
    write.table(myDMP[[1]], file = "myDMP_BAL_HLADR.txt", 
            sep = "\t", 
            quote = FALSE)
}

# Mann-Whitney-Wilcoxon independent sample test
WilcoxonTest <- function(WilcoxonTest){
    B1 <- data1[,c(1,2)]
    combine <- data.frame(c(B1[,"B1_SPU"], B1[, "B1_LAV"]))
    combine$gr <- factor(rep(1:2, each = 398388))
    colnames(combine) <- "combine"
    wilcox.test(combine ~ gr, data = combine)
    wilcox.test(data$B1_SPU_HLDADR, data$B1_LAV_HLDADR, paired = TRUE)
    wilcox.test(combine ~ gr, data = combine)
}

# Hierarchical cluster analysis
ClusterAnalysis <- function(ClusterAnalysis){
    # loading required packages
    suppressMessages(suppressWarnings(library(cluster)))
    suppressMessages(suppressWarnings(library(ape)))

    # Load data from alveolar macrophage (HLA-DR+/CD3-) cells
    hladr <- myNorm[myNorm$HLADR, ]
    HLADR <- t(hladr)
    dist.hladr <- dist(HLADR, method = "euclidean") # calculate the Euclidean distance matrix
    hc.hladr <- hclust(dist.hladr, method = "ward.D2") # cluster analysis using ward method
    colors = c("B1_SPU_HLDADR" = "#FF8C00", "B1_LAV_HLDADR" = "8B0000", 
                "B2_SPU_HLDADR" = "#FF8C00", "B2_LAV_HLDADR" = "8B0000", 
                "B3_SPU_HLDADR" = "#FF8C00", "B3_LAV_HLDADR" = "8B0000", 
                "B4_SPU_HLDADR" = "#FF8C00", "B4_LAV_HLDADR" = "8B0000", 
                "B5_SPU_HLDADR" = "#FF8C00", "B5_LAV_HLDADR" = "8B0000")

    setEPS()
    postscript("hclust_Per_HLADR.eps")
    plot(as.phylo(hc.hladr), 
            font = 1, label.offset = 2,
            no.margin = TRUE)
            tiplabels(pch = 19, col = colors, adj = 0.5, cex = 2)
    dev.off()

    # Load data from alveolar lymphocyte (CD3+/CD3-) cells

    cd3 <- myNorm[myNorm$CD3, ]
    dist.cd3 <- dist(cd3)
    distCD3 <- t(cd3)
    dist.CD3 <- dist(distCD3, method = "euclidean")
    hc_cd3 <- hclust(dist.CD3, method = "ward.D2")
    
    setEPS()
    postscript("hclust_Per_CD3.eps")
    plot(as.phylo(hc_cd3), 
            font = 1, label.offset = 2,
            no.margin = TRUE)
    tiplabels(pch = 19, col = colors, adj = 0.5, cex = 2)
    dev.off()
}