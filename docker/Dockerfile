FROM datashield/rock-base:6.3.0-R4.4.0

USER root

# System dependencies required by units/sf
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    pkg-config \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    && rm -rf /var/lib/apt/lists/*

# units requires Rcpp >= 1.1.0
RUN R -e "install.packages('Rcpp', repos='https://cloud.r-project.org'); stopifnot(packageVersion('Rcpp') >= '1.1.0')"

# Install units first
RUN R -e "install.packages('units', repos='https://cloud.r-project.org'); stopifnot(requireNamespace('units', quietly=TRUE))"

# Install spdep
RUN R -e "install.packages('spdep', repos='https://cloud.r-project.org'); stopifnot(requireNamespace('spdep', quietly=TRUE))"

COPY dsGeospatial_0.0.0.9000.tar.gz /tmp/

RUN R CMD INSTALL /tmp/dsGeospatial_0.0.0.9000.tar.gz

RUN rm /tmp/dsGeospatial_0.0.0.9000.tar.gz
