# The R scripts to analyze the Alveolar macrophages (HLA-DR+/CD3-) and lymphocytes (CD3+) specific cell types from DNA methylation analysis.
[![alv-decon](https://github.com/JD2112/AlveolarCellTypeDeconvolution/actions/workflows/docker-image.yml/badge.svg?event=workflow_run)](https://github.com/JD2112/AlveolarCellTypeDeconvolution/actions/workflows/docker-image.yml)
[![Documentation Status](https://readthedocs.org/projects/alveolarcelltypedeconvolution/badge/?version=latest)](https://alveolarcelltypedeconvolution.readthedocs.io/en/latest/?badge=latest)

## Related publication: (Published in Epigenetics, 2021-08-11)
*Das, J., Idh, N., Paues, J., Sikkeland, L. I. B., & Lerm, M.* (2021). **DNA methylome-based validation of induced sputum as an effective protocol to study lung immunity: construction of a classifier of pulmonary cell types. \ ** bioRxiv.[https://doi.org/10.1101/2021.03.12.435086](https://www.biorxiv.org/content/10.1101/2021.03.12.435086v1) \ [Epigenetics link](https://www.tandfonline.com/doi/full/10.1080/15592294.2021.1969499)

## Create package and R script files according to the analysis (or Result in the manuscript).
1. DNA methylome analysis - till the normalizaed beta value calculation.
2. Normality calculation with Anderson's test (**Table 1**)
3. Pearson's rank correaltion analysis - Figures, Table (**Figure 2 - a. HLA-DR, b. CD3**)
4. Beanplot from the beta values of the whole dataset to describe the beta distribution over all samples (**Figure S1a**)
5. Mann-Whitney test for the hypothesis - Figures, Table (F**igure 3a - HLA-DR and 3b. CD3**)
6. Validation of SI and BAL from Lung compartments (**Figure 4**)
7. Testing of 3 reference-free algorithms - algorithms testings, Venn Diagrams (**Figure 5a. HLA-DR and Figrue 5b. CD3**)
8. Cell proportion analysis using the EpiDISH package (**Figure 6**)

## Use of Docker image
Dockerfile can be used for all R packages and repositories. The image file can be found here 
```
docker pull jd21/alv-decon:latest
```
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
