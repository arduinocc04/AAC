cmake_minimum_required( VERSION 3.8 )
project( no-more-double )

find_library(png
	NAMES libpng
)

find_library(curses
	NAMES ncurses
)

SUBDIRS (
	examples
	external-examples
	img2ascii
)
