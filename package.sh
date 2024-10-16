#!/bin/bash

set -ev

cd /io 

# Build the source distribution
/opt/python/cp38-cp38/bin/python setup.py sdist

echo "Building for Python 3.8"
GRGL_GSL=1 GRGL_BGEN=1 /opt/python/cp38-cp38/bin/python setup.py bdist_wheel
echo "Building for Python 3.9"
GRGL_GSL=1 GRGL_BGEN=1 /opt/python/cp39-cp39/bin/python setup.py bdist_wheel
echo "Building for Python 3.10"
GRGL_GSL=1 GRGL_BGEN=1 /opt/python/cp310-cp310/bin/python setup.py bdist_wheel
echo "Building for Python 3.11"
GRGL_GSL=1 GRGL_BGEN=1 /opt/python/cp311-cp311/bin/python setup.py bdist_wheel
echo "Building for Python 3.12"
GRGL_GSL=1 GRGL_BGEN=1 /opt/python/cp312-cp312/bin/python setup.py bdist_wheel

cd /io/dist
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-1.*-cp38-cp38-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-1.*-cp39-cp39-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-1.*-cp310-cp310-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-1.*-cp311-cp311-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-1.*-cp312-cp312-linux_x86_64.whl
