Alveolar Cell Type Deconvolution
================================

**R scripts to analyze the Alveolar macrophages (HLA-DR+/CD3-) and
lymphocytes (CD3+) specific cell types from DNA methylation analysis.**
.. image:: https://readthedocs.org/projects/alveolarcelltypedeconvolution/badge/?version=latest
    :target: https://alveolarcelltypedeconvolution.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Publication:
------------

| *Das, J., Idh, N., Paues, J., Sikkeland, L. I. B., & Lerm, M.* (2021).
  **DNA methylome-based validation of induced sputum as an effective
  protocol to study lung immunity: construction of a classifier of
  pulmonary cell types.** `Epigenetics link <https://www.tandfonline.com/doi/full/10.1080/15592294.2021.1969499>`__

**BioRxiv.** `https://doi.org/10.1101/2021.03.12.435086 <https://www.biorxiv.org/content/10.1101/2021.03.12.435086v1>`__

Create package and R script files according to the analysis (or Result in the manuscript).
------------------------------------------------------------------------------------------

1. DNA methylome analysis - till the normalizaed beta value calculation.
2. Normality calculation with Anderson’s test (**Table 1**)
3. Pearson’s rank correaltion analysis - Figures, Table (**Figure 2 - a.
   HLA-DR, b. CD3**)
4. Beanplot from the beta values of the whole dataset to describe the
   beta distribution over all samples (**Figure S1a**)
5. Mann-Whitney test for the hypothesis - Figures, Table (F\ **igure 3a
   - HLA-DR and 3b. CD3**)
6. Validation of SI and BAL from Lung compartments (**Figure 4**)
7. Testing of 3 reference-free algorithms - algorithms testings, Venn
   Diagrams (**Figure 5a. HLA-DR and Figrue 5b. CD3**)
8. Cell proportion analysis using the EpiDISH package (**Figure 6**)

Use of Docker image
-------------------

Dockerfile can be used for all R packages and repositories. The image
file can be found here

::

   docker pull jd21/alv-decon:latest

Functions present in the package
--------------------------------

+-------------+------------------+------------------+------------------+
| Functions   | R scripts        | description      | notes            |
+=============+==================+==================+==================+
| *ChAMPanal  | **C              | script for DNA   |                  |
| ysis450K()* | hAMPanalysis.R** | methylation      |                  |
|             |                  | using ChAMP      |                  |
+-------------+------------------+------------------+------------------+
| *Sta        | **Statist        |                  |                  |
| tiscalAnaly | icalAnalysis.R** |                  |                  |
| sisHLADR()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *S          | **Statist        |                  |                  |
| tatiscalAna | icalAnalysis.R** |                  |                  |
| lysisCD3()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *Validatio  | **Vali           |                  |                  |
| nWithCystic | dationWithCF.R** |                  |                  |
| Fibrosis()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *C          | **Statist        |                  |                  |
| ompareAnaly | icalAnalysis.R** |                  |                  |
| sisRingh()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *histo      | **Figure2c.R**   | histogram        |                  |
| gramPlot()* |                  | analysis for     |                  |
|             |                  | beta values      |                  |
+-------------+------------------+------------------+------------------+
| *AlveolarCe | **Figure3.R**    | Houseman         |                  |
| llTypeDecon |                  | algorithm        |                  |
| volutionRef |                  | reference free   |                  |
| FreeEWAS()* |                  | analysis         |                  |
+-------------+------------------+------------------+------------------+
| *Al         | **Figure3.R**    | SVA analysis     |                  |
| veolarCellT |                  |                  |                  |
| ypeDeconvol |                  |                  |                  |
| utionSVA()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *Al         | **Figure3.R**    |                  |                  |
| veolarCellT |                  |                  |                  |
| ypeDeconvol |                  |                  |                  |
| utionRefFre |                  |                  |                  |
| eCellMix()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *Alve       | **Figure3.R**    |                  |                  |
| olarCellTyp |                  |                  |                  |
| eDeconvolut |                  |                  |                  |
| ionTOAST()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *ggplotRe   | **Figure4.R**    |                  |                  |
| gression()* |                  |                  |                  |
+-------------+------------------+------------------+------------------+
| *           | **supplemen      |                  |                  |
| sFigure1()* | taryFigureS1.R** |                  |                  |
+-------------+------------------+------------------+------------------+
| *           | **supplemen      |                  |                  |
| sFigure2()* | taryFigureS2.R** |                  |                  |
+-------------+------------------+------------------+------------------+
| *qqPlot()*  | **supplemen      | Q-Q plot for     | a sub-function   |
|             | taryFigureS3.R** | compare DNA      | can also be      |
|             |                  | methylome data   | used; gg_qq()    |
+-------------+------------------+------------------+------------------+

