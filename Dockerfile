# Example to build:
#   docker build . -t grgl:latest
# Example to run:
#   docker run -v $PWD:/working -it grgl:latest grg_construct /working/test/inputs/msprime.example.vcf

# This file uses a multi-stage docker build. The first stage is the one that grabs the source code and
# builds the GRGL tools and libraries (Python and C++). The results of this build then get copied into
# the second (clean) stage. The second stage is what users will have access to when they run the
# container

# Stage 0: build GRGL
FROM ubuntu:22.04

RUN apt update && \
    apt install -y python3 python3-setuptools python3-pip git build-essential cmake

COPY . /grgl_src

# Install GRGL python API.
RUN cd /grgl_src && pip3 install wheel && GRGL_GSL=1 GRGL_BGEN=1 python3 setup.py bdist_wheel
# Install GRGL command line tools and scripts. Installing the above python package
# will also build this, but there are some extra tools we might want.
RUN cd /grgl_src && mkdir cpp_build && cd cpp_build && mkdir /grgl_inst && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_BGEN=ON -DENABLE_GSL=ON -DCMAKE_INSTALL_PREFIX=/grgl_inst && \
    make -j && \
    make install

# Stage 1: install GRGL and scripts
FROM ubuntu:22.04
COPY --from=0 /grgl_inst/bin/ /usr/local/bin/
COPY --from=0 /grgl_inst/lib/ /usr/local/lib/
COPY --from=0 /grgl_src/dist/pygrgl*.whl /tmp/
COPY --from=0 /grgl_src/scripts/convertlz4.sh /usr/local/bin/gconvertlz4
COPY --from=0 /grgl_src/scripts/experimental/find_treesize.py /usr/local/bin/find_treesize

RUN apt update && \
    apt install -y python3 python3-setuptools python3-pip python-is-python3 lz4 time && \
    pip install --force-reinstall /tmp/pygrgl*.whl && \
    chmod +x /usr/local/bin/gconvertlz4 && \
    chmod +x /usr/local/bin/find_treesize
