#!/usr/bin/env sh

if [ ! -d build ]; then
    echo "build directory doesn't exist. Creating.."
    mkdir build
else
    echo -n "build directory exist. Clear it? If you don't clear build directory, unwanted error may occur. [y/N]"
    read input
    if [[ $input == "" ]] || [[ $input == "n" ]] || [[ $input == "N" ]]; then
        echo "Preserve build directory."
    else
        if [[ $input == "y" ]] || [[ $input == "Y" ]]; then
            rm -rf build/*
            echo "Cleared build directory."
        fi
    fi
fi

cd build

cmake ..

echo make with $(nproc) processes
make -j$(nproc)
echo make finished.

echo generating ascii images
cd ../tools/ascii2image
./install-fonts.sh
