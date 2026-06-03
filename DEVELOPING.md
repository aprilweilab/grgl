# Developing

## Building

Installing via `pip -e` does two things:
1. Any changes to grgl's Python code will be immediately reflected in your current Python environment, without any additional calls to `pip`.
2. When you run the `pip` command below, it compiles the C++ executables and copies them to the top-level of the `grgl` repository. Then when you run the Python code behind `grg construct`, `grg convert`, etc., it will find the executables there. If you make changes to the C++, you'll need to rerun the pip command below; that will rebuild the executables (`grgl`, `grgp`, `grg-merge`) and also rebuild the C++ code behind the Python API.

```
python3 -m venv /path/to/MyEnv
source /path/to/MyEnv/bin/activate
GRGL_COPY_BINS=1 pip install -v -e .
```

Other environment variables that control the Python-based build are:
* `GRGL_DEBUG=1`: Builds the unoptimized version with debug symbols.
* `GRGL_BGEN=1`: Turns on BGEN input support for GRG construction

## Utilities

Use [igdtools](https://pypi.org/project/igdtools/) to perform conversion from `.vcf.gz` to IGD for testing with GRG:
```
pip install igdtools
```

There are also two (less well tested) utilities `gconverter` and `gindexer` that are built when performing the build either via CMake or python. These utilities can be run with `--help` to see options; they are used for converting files (.vcf, .vcf.gz, and BGEN) to IGD and for indexing BGEN files (which is required before providing a BGEN input to GRG construction).

## Testing

### C++ Unit Tests

Build the C++ code manually, in debug mode:
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

Run the tests:
```
./grgl_test  # In build/
```

Add new tests:
* See `test/unit/test_*.cpp` for examples. If you create a new test file, add it to `CMakeLists.txt` in the `GRGL_TEST_SOURCES` list (search for it).

### End-to-end Tests

End-to-end tests can be run via `pytest test/`. These tests build GRGs using the currently installed grgl (so make sure you have built/installed your latest changes).

Adding a new test can be done by adding a method to one of the existing [unittest.TestCase](https://docs.python.org/3/library/unittest.html#unittest.TestCase) classes. That method must be named like `test_*` for pytest to find it.

Adding a new file can be done by creating a file `test/endtoend/test_*.py` (again, so pytest finds it) following a similar pattern to the existing test files.

## Code formatting

* The code is formatted via `clang-format`, and CI workflows on Github will fail unless the code is formatted properly. Run `./format.sh` from the top-level repo directory to format code. `./format-check.sh` can be used to verify the code is properly formatted.
* The `pygrgl/` and `test/endtoend/` directories are formatted via [black](https://pypi.org/project/black/), e.g. `black pygrgl/ test/endtoend/`

