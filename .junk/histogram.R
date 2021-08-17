# Load normalized data
balCD3 <- read.table("myNorm_BAL_CD3.txt", stringsAsFactors = F, header = T)
balHLADR <- read.table("myNorm_BAL_HLADR.txt", stringsAsFactors = F, header = T)

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

#




