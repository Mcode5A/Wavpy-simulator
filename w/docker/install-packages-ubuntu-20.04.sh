#!/bin/sh
apt update && \
apt install -y apt-utils \
               cmake \
               g++ \
               gfortran \
               libgsl23 \
               libgsl-dev \
               libfftw3-bin \
               libfftw3-dev \
               libboost-test-dev \
               doxygen \
               graphviz \
               python3-matplotlib \
               python3-numpy \
               python3-dev && \
adduser --uid 1000 --no-create-home --disabled-password --quiet wavpy

