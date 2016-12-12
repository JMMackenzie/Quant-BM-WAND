#!/bin/bash
cd src
make clean
cd ..
echo "rm -rf build/*"
rm -rf build/*
echo "rm -rf bin/*"
rm -rf bin/*
echo "rm -rf external/*"
rm -rf external/*
