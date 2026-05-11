#!/bin/bash

set -ev

cd /io 

# Build the source distribution
/opt/python/cp39-cp39/bin/pip install setuptools
/opt/python/cp39-cp39/bin/python setup.py sdist

echo "Building for Python 3.9"
/opt/python/cp39-cp39/bin/python setup.py bdist_wheel
echo "Building for Python 3.10"
/opt/python/cp310-cp310/bin/pip install setuptools
/opt/python/cp310-cp310/bin/python setup.py bdist_wheel
echo "Building for Python 3.11"
/opt/python/cp311-cp311/bin/pip install setuptools
/opt/python/cp311-cp311/bin/python setup.py bdist_wheel
echo "Building for Python 3.12"
/opt/python/cp312-cp312/bin/pip install setuptools
/opt/python/cp312-cp312/bin/python setup.py bdist_wheel
echo "Building for Python 3.13"
/opt/python/cp313-cp313/bin/pip install setuptools
/opt/python/cp313-cp313/bin/python setup.py bdist_wheel
echo "Building for Python 3.14"
/opt/python/cp314-cp314/bin/pip install setuptools
/opt/python/cp314-cp314/bin/python setup.py bdist_wheel
echo "Building for Python 3.15"
/opt/python/cp315-cp315/bin/pip install setuptools
/opt/python/cp315-cp315/bin/python setup.py bdist_wheel

cd /io/dist
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp39-cp39-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp310-cp310-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp311-cp311-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp312-cp312-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp313-cp313-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp314-cp314-linux_x86_64.whl
auditwheel repair --plat manylinux_2_24_x86_64 pygrgl-2.*-cp315-cp315-linux_x86_64.whl
