add_executable(hull convexHull.cpp)
target_link_libraries(hull eg)
target_include_directories(hull PRIVATE ../../lib/eg)
