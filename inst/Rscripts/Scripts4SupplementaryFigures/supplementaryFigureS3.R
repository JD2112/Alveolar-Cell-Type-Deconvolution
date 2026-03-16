# Author details
## Script name: qqPlot.R
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


# Q-Q plot from HLA-DR unique CpGs Vs All

qqPlot <- function(qqPlot){
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(ggpubr)))

  ## Loading our data & Stockholm data
  mySIhladr <- merge(uniq_HLADR_BAL, hladr_siAVG, by ="ProbeID")
  stockholmData <- read.table("~/Documents/ML/BALvsIS/Part2- ValidationWithStockholmData/DataSheets/Stockholm_normalizedValue_GSE133062.txt", 
                            stringsAsFactors = FALSE, header = TRUE)
  ## merging CpG data
  st.Data <- merge(uniq_HLADR_BAL, stockholmData, by = "ProbeID")
  colnames(st.Data)[3] <- "uniqHLADR"
  st.Data$stAVG <- rowMeans(st.Data[c(5:39)])

  myData.HLADR <- merge(st.Data, uniq_HLADR_BAL, by ="ProbeID")
  mySI.hladr.stData <- merge(mySIhladr, stockholmData, by = "ProbeID")
  colnames(mySI.hladr.stData)[5] <- "SIavg"
  mySI.hladr.stData$stAVG <- rowMeans(mySI.hladr.stData[c(7,41)])
  # Use non-redundant list of CpGs
  uniq_HLADR_BAL <- merge(uniqCpGsHLADR, hladr_balAVG, by ="ProbeID")

  gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
      q.function <- eval(parse(text = paste0("q", distribution)))
      d.function <- eval(parse(text = paste0("d", distribution)))
      x <- na.omit(x)
      ord <- order(x)
      n <- length(x)
      P <- ppoints(length(x))
      df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
      if(is.null(line.estimate)){
          Q.x <- quantile(df$ord.x, c(0.25, 0.75))
          Q.z <- q.function(c(0.25, 0.75), ...)
          b <- diff(Q.x)/diff(Q.z)
          coef <- c(Q.x[1] - b * Q.z[1], b)
      } else {
          coef <- coef(line.estimate(ord.x ~ z))
        }
  
    zz <- qnorm(1 - (1 - conf)/2)
    SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
    fit.value <- coef[1] + coef[2] * df$z
    df$upper <- fit.value + zz * SE
    df$lower <- fit.value - zz * SE
  
    if(!is.null(labels)){ 
        df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
    }
  
    p <- ggplot(df, aes(x=z, y=ord.x)) +
      geom_point() + 
      geom_abline(intercept = coef[1], slope = coef[2]) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
      theme(text = element_text(family = "Times New Roman", face = "bold", size = 14), 
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 12), axis.line = element_line(size = 1),
            panel.grid.minor = element_blank(),
          )
    if(!is.null(labels)) p <- p + geom_text( aes(label = label))
    print(p)
    coef
  }
  # get the plots
  mod.lm <- lm(log(st.Data$uniqHLADR) ~ log(st.Data$stAVG))
  x <- rstudent(mod.lm)
  gg_qq(x)

  mod.lm1 <- lm(log(mySI.hladr.stData$SIavg) ~ log(mySI.hladr.stData$stAVG))
  y <- rstudent(mod.lm1)
  gg_qq(y)
}