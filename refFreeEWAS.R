# Different Deconvolution method to identify CpGs from heterogeneous samples.
# Rscript modified and written by Jyotirmoy Das.
# Please read the README file for the sample preparation of the analysis.

data <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors=F, header=T)
colnames(data) <- substr(colnames(data), start = 1, stop = 6)

# RefFreeEWAS - Houseman et al, 2014
### Houseman EA, Molitor J, Marsit CJ: Reference-free cell mixture adjustments in analysis of DNA methylation data. Bioinformatics (Oxford, England) 2014, 30(10):1431-1439.
install.packages("RefFreeEWAS")
library("RefFreeEWAS")
edata <- as.matrix(log2(data/(1-data)))
head(edata)
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

## SVA (Surrogate Variable Analysis) based on Singular Value Decomposition (SVD) of residuals
### Leek JT, Storey JD: Capturing heterogeneity in gene expression studies by surrogate variable analysis. PLoS genetics 2007, 3(9):e161.
BiocManager::install(c("sva","MASS","limma"))

library(MASS)
library(sva)
library(limma)
mod1<-model.matrix(~x1+x2+x3+x4)
mod01<-model.matrix(~x2+x3+x4)
svobj1= sva(edata,mod1,mod01,n.sv=NULL,method="two-step")
modSv1 = cbind(mod1,svobj1$sv)
fit1 = lmFit(edata,modSv1)
# method=c("ls","robust")
fite1 = eBayes(fit1)
tab1 = topTable(fite1, coef = "x1",number=length(data[,1]), p.val=0.05,adjust = "fdr")
write.table(tab1,"cpg_SVA_HLA-DR",sep="\t", quote = FALSE)

## RefFreeCellMix
dataY <- as.matrix(data)
cell<-RefFreeCellMix(data,mu0=NULL,K=2,iters=9,Yfinal=NULL,verbose=TRUE)
mod2<-model.matrix(~x1+x2+x3+x4+cell$Omega)
fit2<-lmFit(dataY,mod1,method="robust")
fite2<-eBayes(fit2)
tab2 <- topTable(fite2, coef = "x1",number=length(data[,1]), p.val=0.05,adjust = "fdr")
write.table(tab2,file="RefFreeCellMix_HLA-DR.txt", sep="\t", quote = F)

# Reference-free deconvoluation using RefFreeEWAS, TOAST approach
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TOAST")
library(TOAST)
library(RefFreeEWAS)
library(EpiDISH)

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


# Calculating the wilcox test (two.sided)

comH <- data.frame(c(boxH[,"SI"], boxH[, "BAL"]))
colnames(comH) <- "combine"
#combineB1$gr <- factor(rep(1:2, each = nrow(data)))
comH$grp <- factor(rep(c("SI", "BAL"), each = nrow(boxH)))
wilcox.test(combine ~ grp, data = comH)

#	Wilcoxon rank sum test

#data:  combine by grp
#W = 15, p-value = 0.006841
#alternative hypothesis: true location shift is not equal to 0


comC <- data.frame(c(boxC[,"SI"], boxC[, "BAL"]))
colnames(comC) <- "combine"
comC$grp <- factor(rep(c("SI", "BAL"), each = nrow(boxC)))
wilcox.test(combine ~ grp, data = comC)
#	Wilcoxon rank sum test

#data:  combine by grp
#W = 46, p-value = 0.7959
#alternative hypothesis: true location shift is not equal to 0


> t.test(combine ~ grp, data = comH)

	Welch Two Sample t-test

data:  combine by grp
t = -3.1694, df = 17.116, p-value = 0.005568
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.7026187 -0.1411942
sample estimates:
mean in group BAL  mean in group SI 
        0.2387585         0.6606650 

> t.test(combine ~ grp, data = comC)

	Welch Two Sample t-test

data:  combine by grp
t = -0.33204, df = 17.867, p-value = 0.7437
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4151385  0.3018777
sample estimates:
mean in group BAL  mean in group SI 
        0.4402707         0.4969012 

# Heatmap
library(ComplexHeatmap)
hm <- as.matrix(combine)

colnames(hm) <- c("", "HLA-DR", "", "CD3")
Heatmap(hm, column_title_gp = gpar(fontsize = 12, fontfamily = "Times New Roman"),column_names_rot = 0, row_title_gp = gpar(fontsize = 12, fontfamily = "Times New Roman"), heatmap_legend_param = list(title = "est", color_bar = "discrete"))

























