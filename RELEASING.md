# Releasing new GRGL versions

This outlines the steps for releasing a new GRGL version.

## Version numbering

Bug fixes and minor non-breaking changes do not require a version number bump. Yes, this means that folks who installed via `pip` will not get an updated version unless they `--force-reinstall` or remove and then reinstall. If the bug/change is significant, then bump the minor version.

* Minor version number bump occurs for non-breaking changes.
* Major version number bump occurs for breaking changes. Breaking changes are considered non-backwards-compatible file format changes or breaking changes to the Python API. Breaking changes are allowed to the C++ API without bumping the major version number.

Steps:
* Update `include/grgl/version.h`
* Update `setup.py` to have the new version
* After merging to `main` and ensuring all CI passes (and testing prior to that), run the steps for PyPi below

## Packaging for PyPi

Build the package distributions for PyPi. We build a source dist and then Linux binary distributions for recent Python versions.

```
# Remove dist/ to start fresh
rm -rf dist

# Build the source distribution (for MacOS, mostly)
python setup.py sdist

# Build the Python3.8 and Python3.9 packages on an older Docker image
docker build . -f Dockerfile.pkg1 -t grgl_pkg1:latest
docker run -v $PWD/dist:/dist -it grgl_pkg1:latest bash -c "cp /output/*.whl /dist/"

# Build the Python3.10+ packages on a newer Docker image
docker build . -f Dockerfile.pkg2 -t grgl_pkg2:latest
docker run -v $PWD/dist:/dist -it grgl_pkg2:latest bash -c "cp /output/*.whl /dist/"

# Fix file permissions from Docker
sudo chown ddehaas dist/*.whl
sudo chgrp ddehaas dist/*.whl
```

To upload to PyPi, ensure that your dist directory is clean (only has the packages you _just_ built), then upload via twine:
```
python3 -m twine upload dist/*
```
