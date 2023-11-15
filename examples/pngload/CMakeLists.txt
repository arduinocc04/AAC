add_executable(cvtgray cvtGray.cpp)
target_link_libraries(cvtgray pngload)
target_include_directories(cvtgray PRIVATE ../../lib/PNGLoader)

add_executable(edgedetect edgeDetect.cpp)
target_link_libraries(edgedetect pngload)
target_include_directories(edgedetect PRIVATE ../../lib/PNGLoader)

add_executable(blur-gray blur-gray.cpp)
target_link_libraries(blur-gray pngload)
target_include_directories(blur-gray PRIVATE ../../lib/PNGLoader)
