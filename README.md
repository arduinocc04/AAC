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

# How to use
1. `./install.sh`
2. To use structure-based ascii art, move to build directory and `img2ascii/structure/structure ASCII_FOLDER_NAME INPUT.png SCALE_RATIO`
3. To use structure-based ascii art with simulated annealing optimization , move to build directory and `img2ascii/structure/sa-structure ASCII_FOLDER_NAME INPUT.png`
4. To use tone-based ascii art, move to build directory and `img2ascii/tone/tone INPUT.png`
## Attention
structure-based ascii art with simulated annealing optimization(`main-sa.cpp`) has some debug features. To disable that features, undefine `DEBUG-*`.
## Example
- `cd build && img2ascii/structure/structure non-hangul-images INPUT.png 1`
- `cd build && img2ascii/structure/sa-structure non-hangul-images INPUT.png`  
If you used `install.sh`, `non-hangul-images` will be inside your build directory.
