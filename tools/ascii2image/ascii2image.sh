#!/usr/bin/env sh
#depends: grim
#use kitty in fullscreen.
pos="23,26 14x22"
clear
sleep 1s
for c in `cat $1`
do
    echo -e -n "\r$c"
    sleep 0.01
    echo $pos | grim -g - "$2/${c}.png"
done
