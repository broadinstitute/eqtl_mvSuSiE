FROM rocker/tidyverse:latest

ENV PATH=$PATH:/app \
    GOOGLE_CLOUD_CLI_VERSION="454.0.0" \
    GCLOUD_DIR=/opt/gcloud

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

# Install gsutil with compiled crcmod
RUN apt-get update \
 && apt install curl \
 && mkdir -p $GCLOUD_DIR \
 && curl -so $GCLOUD_DIR/google-cloud-cli.tar.gz \
    https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GOOGLE_CLOUD_CLI_VERSION}-linux-x86_64.tar.gz \
 && tar -xzf $GCLOUD_DIR/google-cloud-cli.tar.gz -C $GCLOUD_DIR \
 && .$GCLOUD_DIR/google-cloud-sdk/install.sh --install-python false --usage-reporting false --rc-path $HOME/.bashrc \
 # manual symlinks since .bashrc is not getting sourced and PATH is getting overwritten in Terra Cloud Environment
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/gsutil /usr/bin/gsutil \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/anthoscli /usr/bin/anthoscli \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/bq /usr/bin/bq \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/docker-credential-gcloud /usr/bin/docker-credential-gcloud \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/gcloud /usr/bin/gcloud

WORKDIR /app