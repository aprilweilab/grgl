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

## Testing

* C++ unit tests are built as `grgl_test`, and running that executable will tell you the status of those tests.
* An end-to-end test can be run via `pytest test/endtoend/run_tests.py`. These tests build GRGs using the currently installed grgl (so make sure you have built/installed your latest changes). The GRGs are then compared against known allele counts to make sure that the mutation to sample relationship was not broken by changes.

## Code formatting

* The code is formatted via `clang-format`, and CI workflows on Github will fail unless the code is formatted properly. Run `./format.sh` from the top-level repo directory to format code. `./format-check.sh` can be used to verify the code is properly formatted.

