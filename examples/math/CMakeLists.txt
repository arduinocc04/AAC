add_executable(hull convexHull.cpp)
target_link_libraries(hull pngload)
target_include_directories(hull PRIVATE ../../lib/PNGLoader)
