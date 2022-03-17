************************************************************************************************************************************************
The R scripts to analyze the Alveolar macrophages (HLA-DR+/CD3-) and lymphocytes (CD3+) specific cell types from DNA methylation analysis.
************************************************************************************************************************************************

Related publication: (Published in Epigenetics, 2021-08-11) 

*Das, J., Idh, N., Paues, J., Sikkeland, L. I. B., & Lerm, M.* (2021). **DNA methylome-based validation of induced sputum as an effective protocol to study lung immunity: construction of a classifier of pulmonary cell types.  **
bioRxiv.`https://doi.org/10.1101/2021.03.12.435086 <https://www.biorxiv.org/content/10.1101/2021.03.12.435086v1>` \
`Epigenetics link <https://www.tandfonline.com/doi/full/10.1080/15592294.2021.1969499>`

Create package and R script files according to the analysis (or Result in the manuscript).
------------------------------------------------------------------------------------------

- DNA methylome analysis - till the normalizaed beta value calculation.
- Normality calculation with Anderson’s test (**Table 1**)
- Pearson’s rank correaltion analysis - Figures, Table (**Figure 2 - a.
   HLA-DR, b. CD3**)
- Beanplot from the beta values of the whole dataset to describe the
   beta distribution over all samples (**Figure S1a**)
- Mann-Whitney test for the hypothesis - Figures, Table (F\ **igure 3a
   - HLA-DR and 3b. CD3**)
- Validation of SI and BAL from Lung compartments (**Figure 4**)
- Testing of 3 reference-free algorithms - algorithms testings, Venn
   Diagrams (**Figure 5a. HLA-DR and Figrue 5b. CD3**)
- Cell proportion analysis using the EpiDISH package (**Figure 6**)

Use of Docker image
-------------------

Dockerfile can be used for all R packages and repositories. The image
file can be found here::

   `docker pull jd21/alv-decon:latest`

Functions present in the package
--------------------------------

+-----------------+-----------------+-----------------+-----------------+
| Functions       | R scripts       | description     | notes           |
+=================+=================+=================+=================+
| *ChAMP          | **Ch            | script for DNA  |                 |
| analysis450K()* | AMPanalysis.R** | methylation     |                 |
|                 |                 | using ChAMP     |                 |
+-----------------+-----------------+-----------------+-----------------+
| *StatiscalA     | **Statisti      |                 |                 |
| nalysisHLADR()* | calAnalysis.R** |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *Statisca       | **Statisti      |                 |                 |
| lAnalysisCD3()* | calAnalysis.R** |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *V              | **Valid         |                 |                 |
| alidationWithCy | ationWithCF.R** |                 |                 |
| sticFibrosis()* |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *CompareA       | **Statisti      |                 |                 |
| nalysisRingh()* | calAnalysis.R** |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *h              | **Figure2c.R**  | histogram       |                 |
| istogramPlot()* |                 | analysis for    |                 |
|                 |                 | beta values     |                 |
+-----------------+-----------------+-----------------+-----------------+
| *AlveolarCellT  | **Figure3.R**   | Houseman        |                 |
| ypeDeconvolutio |                 | algorithm       |                 |
| nRefFreeEWAS()* |                 | reference free  |                 |
|                 |                 | analysis        |                 |
+-----------------+-----------------+-----------------+-----------------+
| *Alveo          | **Figure3.R**   | SVA analysis    |                 |
| larCellTypeDeco |                 |                 |                 |
| nvolutionSVA()* |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *A              | **Figure3.R**   |                 |                 |
| lveolarCellType |                 |                 |                 |
| DeconvolutionRe |                 |                 |                 |
| fFreeCellMix()* |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *Alveola        | **Figure3.R**   |                 |                 |
| rCellTypeDeconv |                 |                 |                 |
| olutionTOAST()* |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *ggpl           | **Figure4.R**   |                 |                 |
| otRegression()* |                 |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *sFigure1()*    | **supplement    |                 |                 |
|                 | aryFigureS1.R** |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *sFigure2()*    | **supplement    |                 |                 |
|                 | aryFigureS2.R** |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| *qqPlot()*      | **supplement    | Q-Q plot for    | a sub-function  |
|                 | aryFigureS3.R** | compare DNA     | can also be     |
|                 |                 | methylome data  | used; gg_qq()   |
+-----------------+-----------------+-----------------+-----------------+

.. figure:: https://github.com/JD2112/AlveolarCellTypeDeconvolution/actions/workflows/docker-image.yml/badge.svg?event=workflow_run
 
.. target:: https://github.com/JD2112/AlveolarCellTypeDeconvolution/actions/workflows/docker-image.yml
