# The R scripts to analyze the Alveolar macrophages (HLA-DR+/CD3-) and lymphocytes (CD3+) specific cell types from DNA methylation analysis.

## Related publication: (Accepted in Epigenetics, Taylor & Francis)
*Das, J., Idh, N., Paues, J., Sikkeland, L. I. B., & Lerm, M.* (2021). **DNA methylome-based validation of induced sputum as an effective protocol to study lung immunity: construction of a classifier of pulmonary cell types.** bioRxiv.[https://doi.org/10.1101/2021.03.12.435086](https://www.biorxiv.org/content/10.1101/2021.03.12.435086v1)

## Create package and R script files according to the analysis (or Result in the manuscript).
1. DNA methylome analysis - till the normalizaed beta value calculation.
2. Normality calculation with Anderson's test (**Table 1**)
3. Spearman's rank correaltion analysis - Figures, Table (**Figure 2 - a. HLA-DR, b. CD3**)
4. Beanplot from the beta values of the whole dataset to describe the beta distribution over all samples (**Figure S1a**)
5. Mann-Whitney test for the hypothesis - Figures, Table (F**igure 3a - HLA-DR and 3b. CD3**)
6. Validation of SI and BAL from Lung compartments (**Figure 4**)
7. Testing of 3 reference-free algorithms - algorithms testings, Venn Diagrams (**Figure 5a. HLA-DR and Figrue 5b. CD3**)
8. Cell proportion analysis using the EpiDISH package (**Figure 6**)


## Functions present in the package
|Functions|R scripts| description|notes|
|---|---|---|---|
|*ChAMPanalysis450K()*|**ChAMPanalysis.R**|script for DNA methylation using ChAMP||
|*StatiscalAnalysisHLADR()*|**StatisticalAnalysis.R**|||
|*StatiscalAnalysisCD3()*|**StatisticalAnalysis.R**|||
|*ValidationWithCysticFibrosis()*|**ValidationWithCF.R**|||
|*CompareAnalysisRingh()*|**StatisticalAnalysis.R**|||
|*histogramPlot()*|**Figure2c.R**|histogram analysis for beta values||
|*AlveolarCellTypeDeconvolutionRefFreeEWAS()*|**Figure3.R**|Houseman algorithm reference free analysis||
|*AlveolarCellTypeDeconvolutionSVA()*|**Figure3.R**|SVA analysis||
|*AlveolarCellTypeDeconvolutionRefFreeCellMix()*|**Figure3.R**|||
|*AlveolarCellTypeDeconvolutionTOAST()*|**Figure3.R**|||
|*ggplotRegression()*|**Figure4.R**|||
|*sFigure1()*|**supplementaryFigureS1.R**|||
|*sFigure2()*|**supplementaryFigureS2.R**|||
|*qqPlot()*|**supplementaryFigureS3.R**|Q-Q plot for compare DNA methylome data|a sub-function can also be used; gg_qq()|

## Figure 1: Correlation between Sputum Induction (SI) and Bronchoalveolar Lavage (BAL)
### Sample analysis - one sample graphical representation
![Image1](https://drive.google.com/uc?export=view&id=1gulQnhXkIp7X3J4XiXAnoMrrevX9ew5r)

### Sample analysis - All samples graphical representation
![Image2](https://drive.google.com/uc?export=view&id=1teYOi7njyPelL8sczUIF-OdoSc1sJE9b)

## Correlation chart
![Image3](https://drive.google.com/uc?export=view&id=1q72_fFlglusYm5HsVTXfb9XAaKp9Tmje)

## Unsupervised hierarchical cluster analysis
![Image4](https://drive.google.com/uc?export=view&id=1rS5ghdUVyHHGtpwLkRvGJq0eKiLdQAic)

