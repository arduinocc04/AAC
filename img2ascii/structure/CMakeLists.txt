add_executable(structure main.cpp)
target_link_libraries(structure eg)
target_include_directories(structure PRIVATE ../../lib/eg)

add_executable(sa-structure main-sa.cpp)
target_link_libraries(sa-structure eg)
target_include_directories(sa-structure PRIVATE ../../lib/eg)
