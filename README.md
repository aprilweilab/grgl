![](https://github.com/aprilweilab/grgl/actions/workflows/cmake-multi-platform.yml/badge.svg)
![](https://readthedocs.org/projects/grgl/badge/?version=latest)

# Genotype Representation Graph Library (GRGL)

GRGL can be used as a library in both C++ and Python. Support is currently limited to Linux and MacOS.
It contains both an API [(see docs)](https://grgl.readthedocs.io/) and a [set of command-line tools](https://github.com/aprilweilab/grgl/blob/main/GettingStarted.md).

## Building (non-Python)

Make sure you clone with `git clone --recursive`!

If you only intend to use GRGL from C++, or use the command-line tools, you can just build it via `CMake`:
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

See below to install the libraries and tools to your system. It is recommended to install it to a custom location (prefix) since removing packages installed via `make install` is a pain otherwise. Example:
```
mkdir /path/to/grgl_installation/
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/grgl_installation/
make -j4
make install
# There should now be bin/, lib/, etc., directories under /path/to/grgl_installation/
```

## Building (Python)

Make sure you clone with `git clone --recursive`!

Requires Python 3.7 or newer to be installed (including development headers). It is recommended that you build/install in a virtual environment.
```
python3 -m venv /path/to/MyEnv
source /path/to/MyEnv/bin/activate
python setup.py bdist_wheel               # Compiles C++, builds a wheel in the dist/ directory
pip install --force-reinstall dist/*.whl  # Install from wheel
```

Or for development you can install the folder with pip:
```
python3 -m venv /path/to/MyEnv
source /path/to/MyEnv/bin/activate
GRGL_COPY_BINS=1 pip install -v -e .
```

BGEN support is disabled by default. If you want to enable it:
* `GRGL_BGEN=1 python setup.py bdist_wheel `
* or `GRGL_BGEN=1 GRGL_COPY_BINS=1 pip install -v -e .`

Build and installation should take at most a few minutes on the typical computer.

## Building (Docker)

We've included a Dockerfile if you want to use GRGL in a container.

Example to build:
```
docker build . -t grgl:latest
```

Example to run, constructing a GRG from an example VCF file:
```
docker run -v $PWD:/working -it grgl:latest bash -c "cd /working && grg construct /working/test/inputs/msprime.example.vcf
```

## Usage (Command line)

There is a command line tool that is mostly for file format conversion and performing common computations on the GRG. For more flexibility, use the Python or C++ APIs.
After building and installing the Python version, run `grg --help` to see all the command options. Some examples are below.

Convert a [tskit](https://tskit.dev/software/tskit.html) tree-sequence into a GRG. This creates `my_arg_data.grg` from `my_arg_data.trees`:
```
grg convert /path/to/my_arg_data.trees my_arg_data.grg
```

Load a GRG and emit some simple statistics about the GRG itself:
```
grg process stats my_arg_data.grg
```

To construct a GRG from a VCF file, use the `grg_from_vcf.py` script (after building grgl):
```
grg construct --parts 20 -j 1 path/to/foo.vcf
```

Construction for small datasets (such as those included as tests in this repository) should be very fast, a few minutes at most. Really large datasets (such as Biobank-scale) can take on the order of a day when using lots of threads (e.g., 70).

## Usage (Python API)

See the provided [jupyter notebooks](https://github.com/aprilweilab/grgl/tree/main/jupyter) and [GettingStarted.md](https://github.com/aprilweilab/grgl/blob/main/GettingStarted.md) for more examples.
