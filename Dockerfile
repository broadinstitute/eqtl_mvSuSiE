FROM python:3.8
COPY src/requirements.txt .

RUN pip3 install -r requirements.txt