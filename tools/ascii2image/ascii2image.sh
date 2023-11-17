#!/usr/bin/env sh
#depends: grim
#use kitty in fullscreen.

if [[ $3 == "" ]]; then
    echo hex or n!! use program correctly.
    echo If you need to capture special character, you may need hex. Otherwise, use n.
    exit
fi
set -f
pos="22,30 11x15" # using slurp will help you get this.
# `clear && echo IIII && echo III` and experiment..
clear
sleep 1s
for c in `cat $1`
do
    if [[ $3 == hex ]]; then
        if [[ $c == \' ]]; then
            echo -n -e "\r'"
            sleep 0.01
            echo $pos | grim -g - "$2/=27.png"
        else
            if [[ $c == \" ]]; then
                echo -n -e "\r\""
                sleep 0.01
                echo $pos | grim -g - "$2/=22.png"
            else
                if [[ $c == \` ]]; then
                    echo -n -e "\r\`"
                    sleep 0.01
                    echo $pos | grim -g - "$2/=60.png"
                else
                    if [[ $c == \\ ]]; then
                        echo -n -e "\r\\"
                        sleep 0.01
                        echo $pos | grim -g - "$2/=5c.png"
                    else
                        if [[ $c == "*" ]]; then
                            echo -n -e "\r*"
                            sleep 0.01
                            echo $pos | grim -g - "$2/=2a.png"
                        else
                            if [[ $c == \# ]]; then
                                echo -n -e "\r#"
                                sleep 0.01
                                echo $pos | grim -g - "$2/=23.png"
                            else
                                echo -e -n "\r$c"
                                sleep 0.01
                                echo $pos | grim -g - "$2/=`echo -n $c | od -A n -t  x1 | sed 's/ *//g'`.png"
                            fi
                        fi
                    fi
                fi
            fi
        fi
    else
        if [[ $3 == n ]]; then
            echo $pos | grim -g - "$2/${c}.png"
        fi
    fi
done
set +f

