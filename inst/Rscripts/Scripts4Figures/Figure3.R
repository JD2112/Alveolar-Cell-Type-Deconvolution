# Reference free EWAS cell proportion analysis
# Different Deconvolution method to identify CpGs from heterogeneous samples.
# Author details
## Script name: Figure 3.R
## Purpose of the script: To visualize the Quantile-Quantile plot from the DNA methylation data mentioned in the Das et al (2021), Epigenetics
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


AlveolarCellTypeDeconvolutionRefFreeEWAS <- function(AlveolarCellTypeDeconvolutionRefFreeEWAS){
        data <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors=F, header=T)
        colnames(data) <- substr(colnames(data), start = 1, stop = 6)

# RefFreeEWAS - Houseman et al, 2014
### Houseman EA, Molitor J, Marsit CJ: Reference-free cell mixture adjustments in analysis 
#of DNA methylation data. Bioinformatics (Oxford, England) 2014, 30(10):1431-1439.

        #install.packages("RefFreeEWAS")
        suppressMessages(suppressWarnings(library("RefFreeEWAS")))
        edata <- as.matrix(log2(data/(1-data)))

        cov1 <- read.table("phenotype.txt", stringsAsFactors=F, header=F)
        cov2 <- read.table("covariates.txt", stringsAsFactors=F, header=F)
        x1 <- as.numeric(cov1[,2])
        x2 <- as.numeric(cov2[,3])
        x3 <- as.numeric(cov2[,5])
        x4 <- as.numeric(cov2[,6])
        tmp1 <- lm(t(edata) ~ x1)
        tmp2 <- cbind(t(coef(tmp1)), t(resid(tmp1)))
        est <- EstDimRMT(tmp2)
        dim3 <- est$dim
        test1 <- RefFreeEwasModel(edata,cbind(1,x1,x2,x3,x4),dim3)
        testBoot1 <- BootRefFreeEwasModel(test1, 50)
        res1 = summary(testBoot1)
        est1 = res1[, 2, 1, 'mean']
        est.sd1 = res1[, 2, 1, 'sd']
        Q21 = (est1/est.sd1)^2
        numCov = 1
        P01 = 2*(pt(-sqrt(Q21),(204-(dim3+6))))
        fdr.pvalue =  p.adjust(P01, method = "fdr", n = length(P01))
                i = which(fdr.pvalue <= 0.05)
                est_ewas1 = est1[i]
        sd_ewas1 = est.sd1[i]
        P0_ewas1 = P01[i]
        fdr_ewas1 = fdr.pvalue[i]
        cpg_ewas1 = cbind(est_ewas1,sd_ewas1,P0_ewas1,fdr_ewas1)
        write.table(cpg_ewas1,file="cpg_RefFreeEWAS_HLA-DR",sep="\t", quote = F)
}

AlveolarCellTypeDeconvolutionSVA <- function(AlveolarCellTypeDeconvolutionSVA){
## SVA (Surrogate Variable Analysis) based on Singular Value Decomposition (SVD) of residuals
### Leek JT, Storey JD: Capturing heterogeneity in gene expression studies by surrogate variable analysis. PLoS genetics 2007, 3(9):e161.

        #BiocManager::install(c("sva","MASS","limma"))

        suppressMessages(suppressWarnings(library(MASS)))
        suppressMessages(suppressWarnings(library(sva)))
        suppressMessages(suppressWarnings(library(limma)))

        mod1<-model.matrix(~x1+x2+x3+x4)
        mod01<-model.matrix(~x2+x3+x4)
        svobj1= sva(edata,mod1,mod01,n.sv=NULL,method="two-step")
        modSv1 = cbind(mod1,svobj1$sv)
        fit1 = lmFit(edata,modSv1)
        # method=c("ls","robust")
        fite1 = eBayes(fit1)
        tab1 = topTable(fite1, coef = "x1",number=length(data[,1]), p.val=0.05,adjust = "fdr")
        write.table(tab1,"cpg_SVA_HLA-DR",sep="\t", quote = FALSE)
}

AlveolarCellTypeDeconvolutionRefFreeCellMix <- function(AlveolarCellTypeDeconvolutionRefFreeCellMix){
        ## RefFreeCellMix
        dataY <- as.matrix(data)
        cell<-RefFreeCellMix(data,mu0=NULL,K=2,iters=9,Yfinal=NULL,verbose=TRUE)
        mod2<-model.matrix(~x1+x2+x3+x4+cell$Omega)
        fit2<-lmFit(dataY,mod1,method="robust")
        fite2<-eBayes(fit2)
        tab2 <- topTable(fite2, coef = "x1",number=length(data[,1]), p.val=0.05,adjust = "fdr")
        write.table(tab2,file="RefFreeCellMix_HLA-DR.txt", sep="\t", quote = F)
}

AlveolarCellTypeDeconvolutionTOAST <- function(AlveolarCellTypeDeconvolutionTOAST){
        # Reference-free deconvoluation using RefFreeEWAS, TOAST approach
        if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")

        BiocManager::install("TOAST")
        suppressMessages(suppressWarnings(library(TOAST)))
        suppressMessages(suppressWarnings(library(RefFreeEWAS)))
        suppressMessages(suppressWarnings(library(EpiDISH)))

        # HLA-DR cell analysis
        Y_h <- as.matrix(dataH)
        refinh <- findRefinx(Y_h, nmarker = 1000)
        Yh <- Y_h[refinh,]
        K <- 2
        outTH <- RefFreeCellMix(Yh, mu0=RefFreeCellMixInitialize(Yh, K = K))
        estProp_RH <- outTH$Omega
        boxH <- estProp_RH

        # CD3 cell analysis
        Y_raw <- as.matrix(dataC)
        refinx <- findRefinx(Y_raw, nmarker = 1000)
        Y <- Y_raw[refinx,]
        K <- 2
        outT <- RefFreeCellMix(Y, mu0=RefFreeCellMixInitialize(Y, K = K))
        estProp_RF <- outT$Omega
        boxC <- estProp_RF

        # Plot them
        combine <- cbind(boxH, boxC)
        boxplot(combine, xaxt='n', frame=FALSE, col = c("#FDD4D4", "#99E5E0", "#FDD4D4", "#99E5E0"))
        axis(side = 1,at = 0:5,labels=c("","HLA-DR", "", "CD3","","" ),lwd.ticks = FALSE)
}

# functions need to be called 
##1. AlveolarCellTypeDeconvolutionRefFreeEWAS() 
##2. AlveolarCellTypeDeconvolutionSVA()
##3. AlveolarCellTypeDeconvolutionRefFreeCellMix()
##4. AlveolarCellTypeDeconvolutionTOAST()