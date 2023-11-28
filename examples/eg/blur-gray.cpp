/**
 * @file blur-gray.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include <iostream>
#include "egLoader.hpp"
#include "egProcessing.hpp"

using namespace eg::imgproc;

int main(int argc, char * argv[]) {
    if(argc != 3) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }
    eg::PNG png;
    png.openImage(argv[1]);
    Image i = png.copy();
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    t = getEdge(t, eg::edgeDetectMethod::gradient);
    int n = std::stoi(argv[2]);
    for(int i = 0; i < n; i++) t = blur(t, eg::blurMethod::gaussian);
    i = mat2dToImage(t);
    png.setImage(i);
    png.saveImage("ed+blur" + std::to_string(n) + "-" + std::string(argv[1]));
}
