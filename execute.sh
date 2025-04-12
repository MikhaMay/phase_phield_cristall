g++-14 \
    -std=c++20 \
    -O3 -march=native -ffast-math -ftree-vectorize \
    code_sources/calculate/main.cpp \
    -o code_sources/calculate/main

./code_sources/calculate/main
./venv/bin/python3.12 code_sources/visualize/__main__.py

