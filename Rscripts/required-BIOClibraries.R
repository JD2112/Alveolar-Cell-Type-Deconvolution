reqd_biocpkgs <- c('ChAMP', 
    'EpiDISH', 
    'limma',
    'TOAST')
requireNamespace("BiocManager")
BiocManager::install(reqd_biocpkgs,ask=F)