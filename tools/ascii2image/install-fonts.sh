#!/usr/bin/env sh
#depends: choose
if [[ `pwd | choose -f / -d -1` != "ascii2image" ]]; then
    echo Use program in ascii2image directory.
    exit
fi

for f in *.tar.gz
do
    fn=`echo -n $f | sed "s/.tar.gz$//"`
    tar xvf $f
    cp -r $fn ../../build
done
