FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgit2-dev

RUN R -e "install.packages('remotes')"

# Install latest versions from CRAN (close to yours)
RUN R -e "install.packages(c('Seurat', 'harmony', 'dplyr', 'ggplot2', 'patchwork'))"
RUN R -e "remotes::install_github('sqjin/CellChat')"


COPY scripts/04_R_analysis/02_integration_harmony.R /app/scripts/
COPY scripts/04_R_analysis/03_cell_type_annotation.R /app/scripts/
COPY scripts/04_R_analysis/04_purkinje_subtype_analysis.R /app/scripts/
COPY scripts/04_R_analysis/05_cellchat_analysis.R /app/scripts/

WORKDIR /app
CMD ["bash"]
