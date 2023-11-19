add_executable(structure main.cpp)
target_link_libraries(structure pngload curses)
target_include_directories(structure PRIVATE ../../lib/PNGLoader)
