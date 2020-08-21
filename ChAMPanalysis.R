# Rscript modified and re-written by Jyotirmoy Das.
# Please read the README file for the sample preparation of the analysis.

library(ChAMP)
myLoad <- champ.load()
myLoad <- champ.load()
myLoad <- champ.load()
myNorm <- champ.norm(beta=myLoad$beta, arraytype = "450K", cores = 30)
myDMP <- champ.DMP(beta = myNorm, pheno=myLoad$pd$Sample_Group)
myDMP <- champ.DMP(beta = myNorm, pheno=myLoad$pd$Sample_Group, adjPVal = 1)
write.table(myNorm, file = "myNorm_BAL_HLADR.txt", sep = "\t", quote = FALSE)
write.table(myDMP[[1]], file = "myDMP_BAL_HLADR.txt", sep = "\t", quote = FALSE)

# Mann-Whitney-Wilcoxon independent sample test
B1 <- data1[,c(1,2)]
combine <- data.frame(c(B1[,"B1_SPU"], B1[, "B1_LAV"]))
combine$gr <- factor(rep(1:2, each = 398388))
colnames(combine) <- "combine"
wilcox.test(combine ~ gr, data = combine)
wilcox.test(data$B1_SPU_HLDADR, data$B1_LAV_HLDADR, paired = TRUE)

wilcox.test(combine ~ gr, data = combine)
