add_library(pngload PNGLoader.cpp PNGMath.cpp edgeDetect.cpp cvtGray.cpp blur.cpp)
set_target_properties(pngload PROPERTIES VERSION 0.0.1)
target_link_libraries(pngload png Eigen3::Eigen)
