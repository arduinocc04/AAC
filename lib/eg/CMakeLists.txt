add_library(eg
    egGeometry.cpp
    egLoader.cpp
    egMath.cpp
    egOperators.cpp
    egProcessing.cpp
    egTrace.cpp
)
set_target_properties(eg PROPERTIES VERSION 0.0.1)
target_link_libraries(eg png Eigen3::Eigen)
