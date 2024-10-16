# Releasing new GRGL versions

This outlines the steps for releasing a new GRGL version.

## Version numbering

Bug fixes and minor non-breaking changes do not require a version number bump. Yes, this means that folks who installed via `pip` will not get an updated version unless they `--force-reinstall` or remove and then reinstall. If the bug/change is significant, then bump the minor version.

* Minor version number bump occurs for non-breaking changes.
* Major version number bump occurs for breaking changes. Breaking changes are considered non-backwards-compatible file format changes or breaking changes to the Python API. Breaking changes are allowed to the C++ API without bumping the major version number.

Steps:
* Update `include/grgl/version.h`
* After merging to `main` and ensuring all CI passes (and testing prior to that), run the steps for PyPi below

## Packaging for PyPi

Build the package distributions for PyPi. We build a source dist and then Linux binary distributions for recent Python versions. The container is based on the [manylinux](https://github.com/pypa/manylinux) project.

```
# Remove dist/ to start fresh
rm -rf dist

# Build the container based on the manylinux project
docker build . -f Dockerfile.package -t grgl_pkg:latest

# Run the packaging inside the container
docker run -v $PWD:/io -it grgl_pkg:latest /io/package.sh

# Fix file permissions from Docker
sudo chown -R ddehaas dist/
sudo chgrp -R ddehaas dist/

# Copy the source wheel to wheelhouse
cp ./dist/*.tar.gz ./dist/wheelhouse/

```

To upload to PyPi:
```
python3 -m twine upload dist/wheelhouse/*
```
