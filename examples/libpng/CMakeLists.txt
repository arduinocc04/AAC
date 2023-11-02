add_executable(readSize readSize.cpp)
target_link_libraries(readSize png)

add_executable(rgb2gray rgb2gray.cpp)
target_link_libraries(rgb2gray png)

add_executable(clearImage clearImage.cpp)
target_link_libraries(clearImage png)
