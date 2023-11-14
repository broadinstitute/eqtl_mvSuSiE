FROM rocker/tidyverse:latest

RUN apt-get update \
    && apt-get install -y libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

ENV MAKE="make -j8"
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("Rfast")'
RUN Rscript -e 'install.packages("mashr")'
RUN Rscript -e 'install.packages("susieR")'
RUN Rscript -e 'remotes::install_github("stephenslab/mvsusieR")'

WORKDIR /app