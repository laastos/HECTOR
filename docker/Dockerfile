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

# Install python development headers for Cython compilation
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy and compile Cython modules
COPY code/sim.pyx /tmp/build/
WORKDIR /tmp/build
RUN cython sim.pyx && \
    gcc -shared -pthread -fPIC -fwrapv -O2 -Wall \
        -fno-strict-aliasing \
        -I/usr/local/include/python3.9 \
        -I$(python3 -c "import numpy; print(numpy.get_include())") \
        -o sim.cpython-39-x86_64-linux-gnu.so \
        sim.c \
        -lpython3.9

# Create entrypoint script to copy compiled modules to workspace
RUN echo '#!/bin/bash\n\
if [ ! -f /workspace/code/sim.cpython-39-x86_64-linux-gnu.so ]; then\n\
    echo "Installing pre-compiled Cython modules..."\n\
    cp /tmp/build/*.so /workspace/code/ 2>/dev/null || true\n\
fi\n\
exec "$@"' > /entrypoint.sh && chmod +x /entrypoint.sh

WORKDIR /workspace

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]
