add_executable(img2ascii-structure main.cpp)
target_link_libraries(img2ascii-structure pngload curses)
target_include_directories(img2ascii-structure PRIVATE ../../lib/PNGLoader)
