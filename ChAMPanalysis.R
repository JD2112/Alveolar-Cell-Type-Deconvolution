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


# Cluster analysis
hladr <- myNorm
head(hladr)
HLADR <- t(hladr)
dist.hladr <- dist(HLADR)
library(cluster)
library(ape)
hc.hladr <- hclust(dist.hladr)
plot(as.phylo(hc.hladr), 
font = 1, label.offset = 20,
no.margin = TRUE)
plot(as.phylo(hc_cd3), 
font = 1, label.offset = 1,
no.margin = TRUE)
head(myNorm, n = 1)
colors = c("B1_SPU_HLDADR" = "#FF8C00", "B1_LAV_HLDADR" = "8B0000", "B2_SPU_HLDADR" = "#FF8C00", "B2_LAV_HLDADR" = "8B0000", "B3_SPU_HLDADR" = "#FF8C00", "B3_LAV_HLDADR" = "8B0000", "B4_SPU_HLDADR" = "#FF8C00", "B4_LAV_HLDADR" = "8B0000", "B5_SPU_HLDADR" = "#FF8C00", "B5_LAV_HLDADR" = "8B0000")
plot(as.phylo(hc.hladr), 
font = 1, label.offset = 1,
no.margin = TRUE)
tiplabels(pch = 19, col = colors, adj = 1, cex = 2)
colors = c("B1_SPU_HLDADR" = "#FF8C00", "B1_LAV_HLDADR" = "#68228B", "B2_SPU_HLDADR" = "#FF8C00", "B2_LAV_HLDADR" = "#68228B", "B3_SPU_HLDADR" = "#FF8C00", "B3_LAV_HLDADR" = "#68228B", "B4_SPU_HLDADR" = "#FF8C00", "B4_LAV_HLDADR" = "#68228B", "B5_SPU_HLDADR" = "#FF8C00", "B5_LAV_HLDADR" = "#68228B")
plot(as.phylo(hc1.hladr), 
font = 1, label.offset = 20,
no.margin = TRUE)
plot(as.phylo(hc.hladr), 
font = 1, label.offset = 1,
no.margin = TRUE)
tiplabels(pch = 19, col = colors, adj = 1, cex = 2)
plot(as.phylo(hc.hladr), 
font = 1, label.offset = 1,
no.margin = TRUE)
tiplabels(pch = 19, col = colors, adj = 2, cex = 2)
plot(as.phylo(hc.hladr), 
font = 1, label.offset = 2,
no.margin = TRUE)
tiplabels(pch = 19, col = colors, adj = 0.5, cex = 2)
setEPS()
postscript("hclust_Per_HLADR.eps")
plot(as.phylo(hc.hladr), 
font = 1, label.offset = 2,
no.margin = TRUE)
tiplabels(pch = 19, col = colors, adj = 0.5, cex = 2)
dev.off()
