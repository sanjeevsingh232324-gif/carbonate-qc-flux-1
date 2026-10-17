FROM python:3.11.13-slim-bookworm@sha256:86adf8dbadc3d6e82ee5dd2c74bec2e1c2467cdad47886280501df722372d2e1

RUN pip install --no-cache-dir numpy==1.26.4 pandas==2.2.2 scipy==1.13.1 PyCO2SYS==1.8.3

WORKDIR /app
COPY data/ /app/
