#!/bin/bash

set -ev

clang-format --dry-run -Werror src/*.cpp
clang-format --dry-run -Werror src/*.h
clang-format --dry-run -Werror src/python/*.cpp
clang-format --dry-run -Werror include/grgl/*.h
