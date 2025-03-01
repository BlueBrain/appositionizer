FROM ubuntu:latest AS builder

RUN apt-get update && apt-get install -y \
 build-essential \
 catch2 \
 cmake \
 git \
 libhdf5-mpich-dev \
 librandom123-dev \
 libtbb-dev \
 libyaml-cpp-dev \
 ninja-build \
 python3 \
 && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive --shallow-submodules --depth=1 --single-branch https://github.com/BlueBrain/morphio.git /tmp/morphio \
 && cmake -B /tmp/morphio/build -S /tmp/morphio -G Ninja -DCMAKE_INSTALL_PREFIX=/opt/appo -DBUILD_BINDINGS:BOOL=OFF -DMorphIO_WERROR:BOOL=OFF -DMORPHIO_TESTS:BOOL=OFF \
 && cmake --build /tmp/morphio/build \
 && cmake --install /tmp/morphio/build \
 && rm -rf /tmp/morphio

RUN git clone --recursive --shallow-submodules --depth=1 --single-branch https://github.com/BlueBrain/libsonata.git /tmp/libsonata \
 && env CXX=mpicxx cmake -B /tmp/libsonata/build -S /tmp/libsonata -G Ninja -DCMAKE_INSTALL_PREFIX=/opt/appo -DSONATA_TESTS:BOOL=OFF -DSONATA_VERSION=1.2.3 -DEXTLIB_FROM_SUBMODULES:BOOL=ON \
 && cmake --build /tmp/libsonata/build \
 && cmake --install /tmp/libsonata/build \
 && rm -rf /tmp/libsonata

RUN mkdir -p /tmp/appo
COPY . /tmp/appo/src
WORKDIR /tmp/appo/src

RUN cmake -B /tmp/appo/build -S . -G Ninja -DCMAKE_INSTALL_PREFIX=/opt/appo -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DEXTLIB_FROM_SUBMODULES:BOOL=ON \
 && cmake --build /tmp/appo/build \
 && cmake --install /tmp/appo/build \
 && (cd /tmp/appo/build  && ctest -VV --test-dir=/tmp/appo/build) \
 && rm -rf /tmp/appo/build

FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
 libhdf5-mpich-103-1t64 \
 libtbb12 \
 libyaml-cpp0.8 \
 mpich \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/appo/bin /opt/appo/bin
RUN mkdir -p /opt/appo/lib
COPY --from=builder /opt/appo/lib/*.so* /opt/appo/lib/

ENV PATH="$PATH:/opt/appo/bin"
