#!/bin/bash

set -ev

clang-format -i src/*.cpp
clang-format -i src/*.h
clang-format -i src/python/*.cpp
clang-format -i include/grgl/*.h
