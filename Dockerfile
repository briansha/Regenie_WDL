# Filename: Dockerfile
# Used for
# - briansha/regenie_r_base:4.1.0

FROM ubuntu:18.04

WORKDIR /cromwell_root/

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl \
      zip \
      git \
      unzip \
      gzip \
      g++ \
      make \
      gfortran \
      zlib1g-dev \
      libgfortran4 \
      liblapacke-dev \
      libopenblas-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      gdebi-core \
      ca-certificates

# R 4.1.0
RUN curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-4.1.0_1_amd64.deb
RUN DEBIAN_FRONTEND=noninteractive gdebi -n r-4.1.0_1_amd64.deb
RUN ln -s /opt/R/4.1.0/bin/R /usr/local/bin/R
RUN ln -s /opt/R/4.1.0/bin/Rscript /usr/local/bin/Rscript

# R packages
RUN R -e "install.packages('data.table', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('qqman', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('tidyverse', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
