hladrData <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors = F, header = T)
head(hladrData)

hladr_si <- hladrData[, c(1,2,4,6,8,10)]
hladr_bal <- hladrData[, c(1,3,5,7,9,11)]

# avarge si and bal data
hladr_si$avg <- (hladr_si$B1_SPU_HLDADR+hladr_si$B2_SPU_HLDADR + hladr_si$B3_SPU_HLDADR + 
                   hladr_si$B4_SPU_HLDADR + hladr_si$B5_SPU_HLDADR)/5
head(hladr_si)
hladr_bal$avg <- (hladr_bal$B1_LAV_HLDADR+hladr_bal$B2_LAV_HLDADR + hladr_bal$B3_LAV_HLDADR + 
                    hladr_bal$B4_LAV_HLDADR + hladr_bal$B5_LAV_HLDADR)/5
head(hladr_bal)

# subgroup only cell specific CpGs
cpgs <- read.table("cellSpecificCpGs.txt", stringsAsFactors = F, header = T)
hladr_cpgs <- merge(cpgs, hladr_bal, by = "ProbeID")
head(hladr_cpgs)
hladrCpGs <- hladr_cpgs[,c(8,2)]
head(hladrCpGs)
colnames(hladrCpGs)[2] <- "type"
# Create a second file for MWW test
hladr_siAVG <- hladr_si[,c(1,7)]
hladr_siAVG$type <- rep("SI")
head(hladr_siAVG)
hladr_balAVG <- hladr_bal[,c(1,7)]
hladr_balAVG$type <- rep("orig")
head(hladr_balAVG)
dim(hladr_balAVG)
allHLADRbal <- hladr_balAVG[,c(2,3)]
head(allHLADRbal)
## Create data for qqplot
qqhladr <- rbind(allHLADRbal, hladrCpGs)
qqPlot(qqhladr$avg)

qqnorm(hladr_balAVG$avg, pch = 1, frame = FALSE)
qqline(hladr_balAVG$avg, col = "steelblue", lwd = 2)


library("car")
qqData <- hladr_cpgs[,c(2, 8)]
qqPlot(qqData$avg)

qqPlot(hladr_cpgs$avg)
qqPlot(hladr_siAVG$avg)
