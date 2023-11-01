#!/usr/bin/env sh
#depend: choose
RED='\033[0;31m'
NOCOLOR='\033[0m'
ln=`ls | grep -F .c | wc -l`
i=1

for f in `ls`;
do
	ext=`echo $f | choose -f '\.' 1`
	fnWithoutExt=`echo $f | choose -f '\.' 0`
	if [[ $ext == c ]]; then
		echo -e "${RED}[$i/$ln]${NOCOLOR}Building ${RED}$f${NOCOLOR}"
		if [[ -e $fnWithoutExt ]]; then
			echo -e "${NOCOLOR}Skip ${RED}$f${NOCOLOR}"
		else
			gcc -lpng $f -o `echo $f | choose -f '\.' 0`
			echo -e "${RED}[$i/$ln]${NOCOLOR}Built ${RED}$f${NOCOLOR}"
		fi
		i=$((i+1))
	fi
done
