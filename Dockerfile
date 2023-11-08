FROM python:3.8

ENV PATH=$PATH:/app

COPY src/requirements.txt .
RUN pip3 install -r requirements.txt

WORKDIR /app