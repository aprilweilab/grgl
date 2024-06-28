# Building

Build all C++ code including documentation.

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
doxygen ../Doxyfile.in
DOC_BUILD_DIR=$PWD sphinx-build -c ../doc/ -b html -Dbreathe_projects.GRGL=$PWD/doc/xml ../doc/ $PWD/doc/sphinx/
```
