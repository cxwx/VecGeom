#!/bin/bash
rm -rf build
mkdir build
cd build
cmake ../ -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DBACKEND=Vc -DVECTOR=avx -DVC_ACCELERATION=ON
make
cd ../