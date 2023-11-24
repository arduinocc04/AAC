add_executable(structure main.cpp)
target_link_libraries(structure eg curses)
target_include_directories(structure PRIVATE ../../lib/eg)
