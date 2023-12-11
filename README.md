# AAC
AAC(Ascii Art Converter, pronounced with sufficient amounts of aggressive) is a program that converts png image to ascii art.  
This project is heavily influenced by structure-based ascii art\[1\].

# SOGANG CSE2035
This program was made for term project of CSE2035 at Sogang University.   
We didn't use image processing library such as OpenCV because using libpng was mandatory.
## no more double
- This is our team name.
- We suffered badly from floating-point error mitigation while solving Hour-line drawing assignments.
- 20231610 조다니엘([arduinocc04](https://github.com/arduinocc04))
  - design program, wrote code
- 20221558 박준영([Park-Joonyoung](https://github.com/Park-Joonyoung))
  - found and fixed memory bug, wrote report, made presentation

# Dependencies
- We used libpng, eigen.  
1. `pacman -S libpng eigen`

# eg
- `eg` is our library. It does - png format image handling, image processing
## Implemented
- monotone chain(convex hull)\[2\]
- grassfire transform\[3\]
- image comparision based on log-polar transform\[4\] and Bhattacharya distance\[5\]. It's influenced by the paper structure-based ascii art\[1\].
- Suzuki's topological structural analysis\[6\]
- basic image tracing
- Bresenham's line algorithm\[7\]
- Cohen–Sutherland algorithm(line clipping)\[8\]

# Documentation
- We used doxygen for documentation.
1. `pacman -S doxygen`
2. `cd doc && doxygen doxy.conf`

# Result
Logo of Arch Linux
![Logo of Arch Linux](presentation/Input.png)
```
                                     W_                                    
                                    J"\                                    
                                   `F `,                                   
                                   J   3                                   
                                  `"   `;                                  
                                  /     `_                                 
                                 `"      \                                 
                                 /       `,                                
                                J"        3                                
                               `"         `;                               
                               J           '_                              
                              `"            \                              
                              /             `,                             
                             -"              ^                             
                            `F               "\                            
                            J_                `,                           
                             ^,                3                           
                           e, `,               "\                          
                          J""^,"^_              `,                         
                         `"   "^.`k              ^                         
                        `/       "~2,            "\                        
                        J          ""h,           `,                       
                       -"                          '_                      
                      `"                            \                      
                      /                             `,                     
                     J"                              `_                    
                    `"                                \                    
                   `F                                 `\                   
                   /                                   `,                  
                  -"                                    ^                  
                 `"                                     "\                 
                _/                 __.,                  `,                
                J"               _~"   "e,                `,               
               ."               ."       ~x                1               
              `"               -"         "\               "\              
              /               `"           `\               `,             
             J"               /             `_               '_            
            ."               `"              \                \            
           `"                J               ^                `\           
           /                 ]               `         `._     `,          
          J"                 F               `,         "'%r-_  1          
         ."                  F                "            "+_"~J\         
        `"                   ;                "              "h_"'         
        /                    \               `"                "k_         
       J"                    ]               `                   ',        
      ."                 _.-~"               `~-,                 "^,      
     `"              _.~""                       ""~._              `,     
     /           _.~""                               ""~,_           ',    
    J"        _.="                                       "^._         \    
   ."      _.r"                                             "^._      `\   
  `"    _.~"                                                   "~,     `,  
  /   _~"                                                         "~,   ', 
 J"`.""                                                             "'._"b 
.L~"                                                                   "~J\
""                                                                       "^
```

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

# References
1. https://www.cse.cuhk.edu.hk/~ttwong/papers/asciiart/asciiart.html
2. https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
3. https://en.wikipedia.org/wiki/Grassfire_transform
4. Computer vision : models, learning, and inference / Simon J. D. Prince.(p. 286-287)
5. https://solanian.github.io/posts/histogram-comparison/
6. https://www.sciencedirect.com/science/article/abs/pii/0734189X85900167
7. https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
8. https://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
