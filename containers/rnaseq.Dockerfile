FROM rocker/r-ver:4.4.3

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    && R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")' \
    && R -e 'BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)' \
    && R -e 'BiocManager::install(c("airway", "DESeq2"), ask = FALSE, update = FALSE)' \
    && Rscript -e 'stopifnot(requireNamespace("airway", quietly = TRUE), requireNamespace("DESeq2", quietly = TRUE))'

WORKDIR /pipeline
