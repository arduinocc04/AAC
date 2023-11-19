#!/usr/bin/env sh
#depends: choose

if [[ `pwd | choose -f / -d -1` != "ascii2image" ]]; then
    echo Use program in ascii2image directory.
    exit
fi

for f in *.tar.gz
do
    fn=`echo -n $f | sed "s/.tar.gz$//"`
    tar xf $f

    if [[ -d ../../build/$fn ]]; then
        echo -n "Previous ascii images of $fn exists. Clear it? [y/N]"
        read input
        if [[ $input == "" ]] || [[ $input == "n" ]] || [[ $input == "N" ]]; then
            echo "Preserve ascii images directory. Skip $fn."
            continue
        else
            if [[ $input == "y" ]] || [[ $input == "Y" ]]; then
                rm -rf ../../build/$fn/*
            fi
        fi
    else
        mkdir ../../build/$fn
    fi

    cp -r $fn/* ../../build/$fn/
done
