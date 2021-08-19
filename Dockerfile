#./Dockerfile
FROM rocker/rstudio

LABEL   authors ="Jyotirmoy Das, Ph.D." \
        description = "Docker image containing all requirement for the Alveolar Cell Type Deconvolution pipeline"

# install the package
RUN git clone https://github.com/JD2112/AlveolarCellTypeDeconvolution.git && \
    cd AlveolarCellTypeDeconvolution

WORKDIR /assests

RUN apt-get update \
        && apt-get install -y --no-install-recommends apt-utils \
        && apt-get install -y --no-install-recommends \
        ## Basic deps
        gdb \
        python3-pip \
        libz-dev \
        liblzma-dev \
        libbz2-dev \
        libpng-dev \
        libgit2-dev \
        ## sys deps from bioc_full
        pkg-config \
        fortran77-compiler \
        byacc \
        automake \
        curl \
        ## new libs
        libglpk-dev \
        ## Databases and other software
        sqlite \
        openmpi-bin \
        mpi-default-bin \
        openmpi-common \
        openmpi-doc \
        tcl8.6-dev \
        tk-dev \
        default-jdk \
        imagemagick \
        tabix \
        ggobi \
        graphviz \
        protobuf-compiler \
        jags \
        ## Additional resources
        xfonts-100dpi \
        xfonts-75dpi \
        biber \
        libsbml5-dev \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

## install R CRAN binary packages
RUN install2.r -e \
testthat

## Install remaining packages from source
COPY ./required-CRANlibraries.R /assests/required-CRANlibraries.R
COPY ./required-BIOClibraries.R /assests/required-BIOClibraries.R
RUN Rscript required-CRANlibraries.R
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libfftw3-dev \
    gcc && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'requireNamespace("BiocManager"); BiocManager::install(ask=F);' \
&& Rscript required-BIOClibraries.R

COPY . /assests

CMD ["Rscript", "alvDcon.R"]