#!/bin/bash

# Определение переменной пути
SOURCE_PATH="code_sources_2d"

# Компиляция
g++-14 \
    -std=c++20 \
    -O3 -march=native -ffast-math -ftree-vectorize \
    -I . \
    -o $SOURCE_PATH/calculate/main \
    $SOURCE_PATH/calculate/main.cpp $SOURCE_PATH/calculate/grid_field.cpp

if [ $? -ne 0 ]; then exit 1; fi

# Расчет
./$SOURCE_PATH/calculate/main

if [ $? -ne 0 ]; then exit 1; fi

# Визуализация
./venv/bin/python3.12 $SOURCE_PATH/visualize/__main__.py

if [ $? -ne 0 ]; then exit 1; fi
