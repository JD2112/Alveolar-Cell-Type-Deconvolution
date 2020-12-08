# New analysis with Cystic Fibrosis data
# Compare the BAL samples with CF data to find out the location
# of the sample collection

# 1. Get the regression analysis
# 2. Calculate the residual values
# 3. Select those CpGs who are at distant from the regression line (here: >\0.3\)
# 4. Use the identified CpGs as the test data to locate the lung position

load("~/Desktop/Link to BackUp_All/BALdata/BALdata/ServerData/HLADR/.RData")
View(mod1)

head(mod1$residuals)
head(mod1$fitted.values)
b5_resd <- cbind(B5_CD3_1$data$ProbeI, mod1$residuals)
head(b5_resd)

colnames(b5_resd) <- c("ProbeID", "residuals")
head(b5_resd)
b5_resd <- cbind(B5_CD3_1$data, mod1$residuals)
head(b5_resd)

colnames(b5_resd)[12] <- "residuals"
b5_resd_posMax <- b5_resd[b5_resd$residuals > 0.3,]
b5_resd_negMax <- b5_resd[b5_resd$residuals < -0.3,]
b5_resd_03 <- rbind(b5_resd_posMax, b5_resd_negMax)
head(b5_resd_03)
dim(b5_resd_03)

# Ready the CF dataset
head(CFdata)
sampleTableCF <- read.table("SampleTable_GSE132547.txt", stringsAsFactors = F, header = T)
head(sampleTableCF)
sampleTableCF <- read.table("SampleTable_GSE132547.txt", stringsAsFactors = F, header = F)
head(sampleTableCF)
colnames(sampleTableCF) <- c("ID", "POSITION")
head(sampleTableCF)
head(CFdata)
colnames(CFdata)[2:25] <- sampleTableCF$POSITION
head(CFdata)

newcols <- sapply(c("RLL$", "RUL$"), function(x) rowMeans(CFdata[grep(x, names(CFdata))]))
newcols1 <- setNames(cbind(CFdata[1], newcols), c(names(CFdata)[1], "RLL", "RUL"))

# Perform the analysis
library(EpiDISH)
referenceData <- newcols1
names(referenceData) <- referenceData[,1]
names(referenceData) <- newcols1[,1]
head(referenceData)
rownames(referenceData) <- referenceData[,1]
referenceData <- referenceData[,2:3]
head(referenceData)
refData.m <- as.matrix(referenceData)
methData <- b5_resd_03
head(methData)
rownames(methData) <- methData[,1]
methData <- methData[, 2:11]
methData.m <- as.matrix(methData)
out.l5 <- epidish(beta.m = methData.m, ref.m = refData.m, method = "RPC")$estF
boxplot(out.l5, col = c("#FDD4D4", "#99E5E0"))

# Result from CFdata vs B5 meth data CD3 file
> out.l5
            RLL      RUL
B5_SPU 0.520214 0.479786
B5_LAV 0.000000 1.000000

# Result from CFdata vs B1 meth data HLADR file
> out.l1
                    RLL       RUL
B1_SPU_HLDADR 1.0000000 0.0000000
B1_LAV_HLDADR 0.7902722 0.2097278

# Result from CFdata vs B2 meth data HLADR file
> out.l2
              RLL RUL
B2_SPU_HLDADR   1   0
B2_LAV_HLDADR   1   0

# Result from CFdata vs B3 meth data HLADR file
> out.l3
              RLL RUL
B3_SPU_HLDADR   1   0
B3_LAV_HLDADR   1   0

# Result from CFdata vs B4 meth data HLADR file
> out.l4
              RLL RUL
B4_SPU_HLDADR   0   1
B4_LAV_HLDADR   1   0

# Result from CFdata vs B5 meth data HLADR file
> out.l5
              RLL RUL
B5_SPU_HLDADR   0   1
B5_LAV_HLDADR   1   0

# Result from CFdata vs B1 meth data CD3 file
> out.lcd3b1
           RLL RUL
B1_SPU_CD3   1   0
B1_LAV_CD3   1   0

> out.lcd3b2
           RLL RUL
B2_SPU_CD3   1   0
B2_LAV_CD3   1   0

> out.lcd3b3
           RLL RUL
B2_SPU_CD3   1   0
B2_LAV_CD3   1   0

> out.lcd3b4
           RLL RUL
B4_SPU_CD3   0   1
B4_LAV_CD3   1   0

> out.lcd3b5
                RLL      RUL
B5_SPU_CD3 0.520214 0.479786
B5_LAV_CD3 0.000000 1.000000