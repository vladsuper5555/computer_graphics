#!/bin/bash

# Compile the C++ program with OpenGL libraries
echo "Compiling cg0.cpp..."
g++ -o cg0.bin cg0.cpp -lGL -lGLU -lglut

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful! Running the program..."
    # Execute the compiled binary
    ./cg0.bin
else
    echo "Compilation failed!"
    exit 1
fi
