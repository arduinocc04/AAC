add_executable(tone main.cpp)
target_link_libraries(tone eg)
target_include_directories(tone PRIVATE ../../lib/eg)
