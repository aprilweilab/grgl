# Getting Started

This guide walks you through using GRGL locally on a Linux or MacOS system. Requirements:
* [CMake](https://cmake.org/download/) v3.14 or newer
* GCC or Clang that supports C++11
* Python 3.7 or newer
* zlib/libz development libraries

## Build and Install

**Create a Python virtual environment to work in.**

```
python3 -m venv /path/to/anEnv
```

Replace `/path/to/anEnv` with any directory where you want to create the virtual environment.

**Enable the Python virtual environment**

```
source /path/to/anEnv/bin/active
```

All Python packages you install will now be specific to this environment. You can deactivate the environment
with the command `deactivate`. You can delete the environment directory when you are done with it.

**Build and Install GRGL**

```
git clone --recursive https://github.com/aprilweilab/grgl.git
cd grgl
pip install wheel
python setup.py bdist_wheel --no-bgen
pip install --force-reinstall dist/*.whl
```

## Running an Example

Construct a GRG for one of the test inputs. From within the virtual environment, run:
```
grg construct test/endtoend/input/test-200-samples.vcf.gz
```

This will create the GRG `test-200-samples.vcf.gz.final.grg` in the directory you ran the command from.
Using the Python API, you can load the GRG and calculate the allele counts using:
```
import pygrgl
from pygrgl import get_dfs_order, TraversalDirection

g = pygrgl.load_immutable_grg("test-200-samples.vcf.gz.final.grg")

node_sample_counts = {}
for node_id in get_dfs_order(g, TraversalDirection.DOWN, g.get_root_nodes()):
    count = 1 if g.is_sample(node_id) else 0
    for child_id in g.get_down_edges(node_id):
        count += node_sample_counts[child_id]
    node_sample_counts[node_id] = count

for node_id, mutation_id in g.get_node_mutation_pairs():
    mutation = g.get_mutation_by_id(mutation_id)
    if node_id != pygrgl.INVALID_NODE:
        count = node_sample_counts[node_id]
        print(f"Mutation({mutation_id}): ({mutation.position}, {mutation.allele}) has sample count {count}")
```

This code is also in the included Jupyter Notebook in `grgl/jupyter/01.Intro.ipynb`. To run that notebook, do the following in
your `grgl` directory, in the virtual environment:
```
pip install jupyterlab # Install jupyter lab in your virtual environment if you haven't yet
jupyter lab .
```

In your browser, open the URL that JupyterLab prints to the screen. Navigate the tree to `jupyter/01.Intro.ipynb` and run the
cells in the notebook one at a time.

## Next steps

When working with anything but the smallest datasets, you'll want to convert your data to
[IGD](https://github.com/aprilweilab/picovcf/blob/main/README.md#indexable-genotype-data-igd)
before generating a GRG. This can be done with the `grg convert` command, which can convert
any of `.vcf`, `.vcf.gz`, or `.bgen` to `.igd`.

### Converting data to IGD

```
grg convert <.vcf/.vcf.gz/.bgen file> <.igd output filename>
```

### Converting Tree-Sequence to GRG

```
grg convert <.trees input file> <.grg output filename>
```

### Processing GRG file

1. Allele frequency

```
grg process freq <.grg input filename>
```

2. Graph statistics
```
grg process stats <.grg input filename>
```

