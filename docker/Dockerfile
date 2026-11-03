# ---------------------------------------------------------------------------
# GLOBAL VARIABLES
# Override at build time, e.g.:
#   docker build --build-arg DSGEOSPATIAL_REF=dev-MoransI -t ds-geospatial-rock:dev-MoransI .
#   docker build --build-arg DSGEOSPATIAL_REF=93de9cd     -t ds-geospatial-rock:93de9cd .
# ---------------------------------------------------------------------------
ARG ROCK_BASE_IMAGE=datashield/rock-base:6.3.0-R4.4.0
ARG DSGEOSPATIAL_REPO=FederatedMethods/dsGeospatial
# Branch, tag or commit SHA of the GitHub repo to install
ARG DSGEOSPATIAL_REF=main

FROM ${ROCK_BASE_IMAGE}

# ARGs declared before FROM must be re-declared to be used after it
ARG DSGEOSPATIAL_REPO
ARG DSGEOSPATIAL_REF
# Change this value (e.g. --build-arg CACHEBUST=$(date +%s)) to force Docker to
# re-fetch the latest commit of a moving ref such as a branch
ARG CACHEBUST=0

USER root

# System dependencies required by units/sf
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    pkg-config \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# units requires Rcpp >= 1.1.0
RUN R -e "install.packages('Rcpp', repos='https://cloud.r-project.org'); stopifnot(packageVersion('Rcpp') >= '1.1.0')"

# Install units first
RUN R -e "install.packages('units', repos='https://cloud.r-project.org'); stopifnot(requireNamespace('units', quietly=TRUE))"

# Install spdep
RUN R -e "install.packages('spdep', repos='https://cloud.r-project.org'); stopifnot(requireNamespace('spdep', quietly=TRUE))"

# Install remotes (used to install from GitHub)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org'); stopifnot(requireNamespace('remotes', quietly=TRUE))"

# Install {dsGeospatial} from GitHub at the requested branch/tag/commit.
# CACHEBUST is echoed so that changing it invalidates this layer.
# Dependencies are already installed above, so upgrade = "never" keeps the
# versions pinned by the base image / steps above.
RUN echo "Installing ${DSGEOSPATIAL_REPO}@${DSGEOSPATIAL_REF} (cachebust=${CACHEBUST})" \
    && R -e "remotes::install_github('${DSGEOSPATIAL_REPO}', ref='${DSGEOSPATIAL_REF}', upgrade='never', dependencies=TRUE, force=TRUE); stopifnot(requireNamespace('dsGeospatial', quietly=TRUE)); cat('Installed dsGeospatial', as.character(packageDescription('dsGeospatial')[['GithubSHA1']]), '\n')"

# Record what was installed in the image metadata
LABEL org.opencontainers.image.source="https://github.com/${DSGEOSPATIAL_REPO}" \
      dsgeospatial.ref="${DSGEOSPATIAL_REF}"
