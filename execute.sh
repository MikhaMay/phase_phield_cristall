#!/bin/bash

# Компиляция 
g++-14 \
    -std=c++20 \
    -O3 -march=native -ffast-math -ftree-vectorize \
    -I . \
    -o code_sources/calculate/main \
    code_sources/calculate/main.cpp code_sources/calculate/grid_field.cpp

if [ $? -ne 0 ]; then exit 1; fi


# Расчет
./code_sources/calculate/main

if [ $? -ne 0 ]; then exit 1; fi

# Визуализация
./venv/bin/python3.12 code_sources/visualize/__main__.py

if [ $? -ne 0 ]; then exit 1; fi
