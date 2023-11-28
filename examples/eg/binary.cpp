/**
 * @file binary.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include <iostream>
#include "egLoader.hpp"
#include "egProcessing.hpp"

using namespace eg;
using namespace eg::imgproc;

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }
    PNG png;
    png.openImage(argv[1]);
    Image i = png.copy();
    Mat2d t = cvtGray(i, grayCvtMethod::mean);
    t = getEdge(t, edgeDetectMethod::gradient);
    for(int i = 0; i < 10; i++) {
        std::cout << "\r" << i + 1 << "/10 blurring input image";
        t = blur(t, blurMethod::gaussian);
    }
    t = binary(t, 10);
    i = mat2dToImage(t);
    png.setImage(i);
    png.saveImage("bin-" + std::string(argv[1]));
}
