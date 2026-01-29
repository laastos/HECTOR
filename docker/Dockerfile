FROM python:3.9-slim-bookworm

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        make \
        wget \
        unzip \
    && rm -rf /var/lib/apt/lists/*

# Install EDTSurf for surface mesh generation
RUN wget https://zhanggroup.org/EDTSurf/EDTSurf.zip \
    && unzip EDTSurf.zip \
    && g++ -O3 -o EDTSurf EDTSurf.cpp ParsePDB.cpp ProteinSurface.cpp \
    && mv EDTSurf /usr/local/bin/ \
    && rm -rf EDTSurf.zip *.cpp *.h Makefile

RUN pip install --no-cache-dir \
        cython==3.0.11 \
        joblib==1.4.2 \
        jupyter==1.0.0 \
        matplotlib==3.6.2 \
        numpy==1.26.0 \
        scipy==1.12.0

WORKDIR /workspace

CMD ["/bin/bash"]
