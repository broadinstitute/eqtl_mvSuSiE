FROM python:3.8

ENV PATH=$PATH:/app \
    GOOGLE_CLOUD_CLI_VERSION="454.0.0" \
    GCLOUD_DIR=/opt/gcloud

COPY src/requirements.txt .
RUN pip3 install -r requirements.txt

# Install gsutil with compiled crcmod
RUN mkdir -p $GCLOUD_DIR \
 && curl -so $GCLOUD_DIR/google-cloud-cli.tar.gz \
    https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${GOOGLE_CLOUD_CLI_VERSION}-linux-x86_64.tar.gz \
 && tar -xzf $GCLOUD_DIR/google-cloud-cli.tar.gz -C $GCLOUD_DIR \
 && .$GCLOUD_DIR/google-cloud-sdk/install.sh --usage-reporting false --rc-path $HOME/.bashrc \
 # manual symlinks since .bashrc is not getting sourced and PATH is getting overwritten in Terra Cloud Environment
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/gsutil /usr/bin/gsutil \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/anthoscli /usr/bin/anthoscli \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/bq /usr/bin/bq \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/docker-credential-gcloud /usr/bin/docker-credential-gcloud \
 && ln -s $GCLOUD_DIR/google-cloud-sdk/bin/gcloud /usr/bin/gcloud

WORKDIR /app