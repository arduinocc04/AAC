#!/usr/bin/env sh

find_cmake () {
    for f in `ls $1`
    do
        if [ -d "$1/$f" ]; then
            find_cmake "$1/$f" "$2"
        else
            if [[ $f == "CMakeLists.txt" ]]; then
                if [ ! -d "$2/$1" ]; then
                    mkdir $2/$1
                fi
                cp "$1/$f" "$2/$1"
            fi
        fi
    done
}

if [ ! -d $2 ]; then
    mkdir $2
fi

find_cmake . $1
