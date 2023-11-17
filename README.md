# AAC
AAC(Ascii Art Converter, pronounced with sufficient amounts of aggressive) is a program that converts png image to ascii art.

# SOGANG CSE2035
This program was made for term project of CSE2035 at Sogang University.   
We didn't use image processing library such as OpenCV because using libpng was mandatory.
## no more double
- This is our team name.
- We suffered badly from floating-point error mitigation while solving Hour-line drawing assignments.
- 20231610 조다니엘([arduinocc04](https://github.com/arduinocc04))
- 20221558 박준영([Park-Joonyoung](https://github.com/Park-Joonyoung))

# Dependencies
- We use libpng, ncurses, eigen.  
1. `pacman -S libpng ncurses eigen`

# Documentation
- We use doxygen for documentation.  
1. `pacman -S doxygen`
2. `cd doc && doxygen doxy.conf`

# How to build
1. `mkdir build && cd build`
3. `cmake ..`
4. `make`

# How to use
1. Generate ascii image file: generate ascii character file and `tools/ascii2image/ascii2image.sh ASCII_CHARACTER_FILE OUTPUT_FOLDER`

