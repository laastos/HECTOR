# hash:sha256:72c9f6200e1526cae1372b4402e65301144ee3fba9ed39ba9848df8861d98704
FROM registry.codeocean.com/codeocean/miniconda3:4.12.0-python3.9-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gcc=4:9.3.0-1ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        cython==3.0.11 \
        joblib==1.4.2 \
        jupyter==1.0.0 \
        matplotlib==3.6.2 \
        numpy==1.23.4 \
        scipy==1.12.0 \
    && conda clean -ya
