add_executable(png-load cvtGray.cpp)
target_link_libraries(png-load pngload)
target_include_directories(png-load PRIVATE ../../lib/PNGLoader)
