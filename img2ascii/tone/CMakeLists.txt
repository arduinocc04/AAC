add_executable(tone main.cpp)
target_link_libraries(tone pngload)
target_include_directories(tone PRIVATE ../../lib/PNGLoader)
