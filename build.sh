#!/usr/bin/env sh
set -eu

compiler="${CXX:-g++}"
output="hicreate"

"$compiler" -std=c++17 -O2 -Wall -Wextra -pedantic \
  -o "$output" \
  main.cpp matrix.cpp reference.cpp fragmenter.cpp simulator.cpp
