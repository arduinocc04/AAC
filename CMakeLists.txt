cmake_minimum_required( VERSION 3.8 )
project( no-more-double )

find_library(png
	NAMES libpng
)

SUBDIRS (
	libpng-demo
)
