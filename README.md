![](https://github.com/aprilweilab/grgl/actions/workflows/cmake-multi-platform.yml/badge.svg)
![](https://readthedocs.org/projects/grgl/badge/?version=latest)

# Genotype Representation Graphs

A Genotype Representation Graph (GRG) is a compact way to store reference-aligned genotype data for large
genetic datasets. These datasets are typically stored in tabular formats (VCF, BCF, BGEN, etc.) and then
compressed using off-the-shelf compression. In contrast, a GRG contains Mutation nodes (representing variants)
and Sample nodes (representing haploid samples), where there is a path from a Mutation node to a Sample
node if-and-only-if that sample contains that mutation. These paths go through internal nodes that represent
common ancestry between multiple samples, and this can result in significant compression **(30-50x smaller than
.vcf.gz)**. Calculations on the whole dataset can be performed very quickly on GRG, using GRGL. See our paper
["Enabling efficient analysis of biobank-scale data with genotype representation graphs"](https://www.nature.com/articles/s43588-024-00739-9)
for more details.

Since the publication of the paper, [version 2.0](https://github.com/aprilweilab/grgl/releases/tag/v2.0) has been released,
which further reduced the GRG size (by about half) and significantly sped up graph load time (by about 20x).

# Genotype Representation Graph Library (GRGL)

GRGL can be used as a library in both C++ and Python. Support is currently limited to Linux and MacOS.
It contains both an API [(see docs)](https://grgl.readthedocs.io/) and a [set of command-line tools](https://github.com/aprilweilab/grgl/blob/main/GettingStarted.md).

## Installing from pip

If you just want to use the tools (e.g., constructing GRG or converting tree-sequence to GRG) and the Python API then you can install via pip (from [PyPi](http://pypi.org/project/pygrgl/)).

```
pip install pygrgl
```

This will use prebuilt packages for most modern Linux situations, and will build from source for MacOS. In order to build from source it will require CMake (at least v3.14), zlib development headers, and a clang or GCC compiler that supports C++11.

## Building (Python)

The Python installation installs the command line tools and Python libraries (the C++ executables are packaged as part of this). Make sure you clone with `git clone --recursive`!

Requires Python 3.7 or newer to be installed (including development headers). It is recommended that you build/install in a virtual environment.
```
python3 -m venv /path/to/MyEnv
source /path/to/MyEnv/bin/activate
python setup.py bdist_wheel               # Compiles C++, builds a wheel in the dist/ directory
pip install --force-reinstall dist/*.whl  # Install from wheel
```

Build and installation should take at most a few minutes on the typical computer. For more details on build options, see [DEVELOPING.md](https://github.com/aprilweilab/grgl/blob/main/DEVELOPING.md).

## Building (C++ only)

The C++ build is only necessary for folks who want to include GRGL as a library in their C++ project. Typically, you would include our
CMake into your project via [add\_subdirectory](https://cmake.org/cmake/help/latest/command/add_subdirectory.html), but you can also build
standalone as below. Make sure you clone with `git clone --recursive`!

If you only intend to use GRGL from C++, you can just build it via `CMake`:
```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

See below to install the libraries to your system. It is recommended to install it to a custom location (prefix) since removing packages installed via `make install` is a pain otherwise. Example:
```
mkdir /path/to/grgl_installation/
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/grgl_installation/
make -j4
make install
# There should now be bin/, lib/, etc., directories under /path/to/grgl_installation/
```

## Building (Docker)

We've included a Dockerfile if you want to use GRGL in a container.

Example to build:
```
docker build . -t grgl:latest
```

Example to run, constructing a GRG from an example VCF file:
```
docker run -v $PWD:/working -it grgl:latest bash -c "cd /working && grg construct /working/test/inputs/msprime.example.vcf"
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

To construct a GRG from a VCF file, use the `grg construct` command:
```
grg construct --parts 20 -j 1 path/to/foo.vcf
```

**WARNING:** VCF access for GRG is not indexed, and in general really slow. For anything beyond toy datasets, it is recommended to convert
VCF files to [IGD](https://github.com/aprilweilab/picovcf) first. You can use the `grg convert` tool (available as part of GRGL)
 or `igdtools` from [picovcf](https://github.com/aprilweilab/picovcf).

To convert a VCF(.gz) to an IGD and then build a GRG:
```
grg convert path/to/foo.vcf foo.igd
grg construct --parts 20 -j 1 foo.igd
```

Construction for small datasets (such as those included as tests in this repository) should be very fast, a few minutes at most. Really large datasets (such as Biobank-scale) can take on the order of a day when using lots of threads (e.g., 70).

## Usage (Python API)

See the provided [jupyter notebooks](https://github.com/aprilweilab/grgl/tree/main/jupyter) and [GettingStarted.md](https://github.com/aprilweilab/grgl/blob/main/GettingStarted.md) for more examples.
