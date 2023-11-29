/**
 * @file edgeDetect.cpp
 * @author Daniel Cho
 * @date 23.11.29
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
    t = markOutlier(t, 10);
    int n = std::stoi(argv[2]);
    for(int i = 0; i < n; i++) {
        t = dilate(t, 1, 3);
        t = erode(t, 1, 3);
    }
    for(int i = 0; i < n; i++) {
        t = dilate(t, 3, 1);
        t = erode(t, 3, 1);
    }
    t = 255*t;
    i = mat2dToImage(t);
    png.setImage(i);
    std::string nn = argv[2];
    png.saveImage("morph+" + nn + "-" + std::string(argv[1]));
}
