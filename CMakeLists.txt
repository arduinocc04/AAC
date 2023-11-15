cmake_minimum_required( VERSION 3.8 )
project( no-more-double )

set (CMAKE_CXX_STANDARD 17)

find_library(png
	NAMES libpng
)

find_library(curses
	NAMES ncurses
)

find_package(Eigen3 3.3 REQUIRED)

SUBDIRS (
    lib
	examples
	external-examples
	img2ascii
)
