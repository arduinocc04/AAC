add_executable(cvtgray cvtGray.cpp)
target_link_libraries(cvtgray eg)
target_include_directories(cvtgray PRIVATE ../../lib/eg)

add_executable(edgedetect edgeDetect.cpp)
target_link_libraries(edgedetect eg)
target_include_directories(edgedetect PRIVATE ../../lib/eg)

add_executable(blur-gray blur-gray.cpp)
target_link_libraries(blur-gray eg)
target_include_directories(blur-gray PRIVATE ../../lib/eg)

add_executable(choose-most-special choose-most-special.cpp)
target_link_libraries(choose-most-special eg)
target_include_directories(choose-most-special PRIVATE ../../lib/eg)

add_executable(binary binary.cpp)
target_link_libraries(binary eg)
target_include_directories(binary PRIVATE ../../lib/eg)

add_executable(masking masking.cpp)
target_link_libraries(masking eg)
target_include_directories(masking PRIVATE ../../lib/eg)

add_executable(exCenter extractCenterline.cpp)
target_link_libraries(exCenter eg)
target_include_directories(exCenter PRIVATE ../../lib/eg)
