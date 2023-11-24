add_library(eg
    egLoader.cpp
    egMath.cpp
    egProcessing.cpp
    egGeometry.cpp
    trace.cpp
)
set_target_properties(eg PROPERTIES VERSION 0.0.1)
target_link_libraries(eg png Eigen3::Eigen)
