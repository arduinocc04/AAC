#!/usr/bin/env sh

echo "" > /tmp/count_line_data
exts=( cpp hpp c h sh py md )

count_line () {
    for f in `ls $1`
    do
        if [ -d "$1/$f" ]; then
            if [[ $f != "build" ]] && [[ $f != "debug" ]]; then
                count_line "$1/$f"
            fi
        else
            for ext in "${exts[@]}"
            do
                if [[ $f == *.$ext ]] || [[ $f == "CMakeLists.txt" ]]; then
                    t=`cat $1/$f | wc -l`
                    echo Counted $1/$f $t
                    echo $t >> /tmp/count_line_data
                    break
                fi
            done
        fi
    done
}

count_line .

s=0
for n in `cat /tmp/count_line_data`
do
    s=$((s+n))
done

echo $s
