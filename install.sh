#!/usr/bin/env sh
if [ ! -d build ]; then
    mkdir build
fi

cd build
cmake ..
make -j$(nproc)

cd ../tools/ascii2image
./install-fonts.sh
