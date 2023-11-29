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

add_executable(get-contours get-contours.cpp)
target_link_libraries(get-contours eg)
target_include_directories(get-contours PRIVATE ../../lib/eg)

add_executable(histogram histogram.cpp)
target_link_libraries(histogram eg)
target_include_directories(histogram PRIVATE ../../lib/eg)

add_executable(cycle-image-and-comp cycle-image-and-comp.cpp)
target_link_libraries(cycle-image-and-comp eg)
target_include_directories(cycle-image-and-comp PRIVATE ../../lib/eg)

add_executable(approx-border approx-border.cpp)
target_link_libraries(approx-border eg)
target_include_directories(approx-border PRIVATE ../../lib/eg)

add_executable(morph-gray morph-gray.cpp)
target_link_libraries(morph-gray eg)
target_include_directories(morph-gray PRIVATE ../../lib/eg)

