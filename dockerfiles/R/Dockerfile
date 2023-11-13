FROM rocker/tidyverse:latest

RUN apt-get update

RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("Rfast")'
RUN Rscript -e 'install.packages("mashr")'
RUN Rscript -e 'install.packages("susieR")'
RUN Rscript -e 'remotes::install_github("stephenslab/mvsusieR")'

WORKDIR /app