reqd_biocpkgs <- c('ChAMP', 
    'ggplot2', 
    'ggpubr', 
    'EpiDISH', 
    'limma',
    'TOAST'
)

requiredNamespace("BiocManager")
BiocManager::install(reqd_biocpkgs, ask = F)