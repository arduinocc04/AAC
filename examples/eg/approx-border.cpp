/**
 * @file approx-border.cpp
 * @author Daniel Cho
 * @date 2023.11.29
 * @version 0.0.1
 */
#include <iostream>

#include "egProcessing.hpp"
#include "egTrace.hpp"
#include "egLoader.hpp"

using namespace eg::imgproc;

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use Program properly!" << std::endl;
        return -1;
    }
    std::string inputPath = argv[1];

    PNG png;
    png.openImage(inputPath);
    Image i = png.copy();
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    t = getEdge(t, edgeDetectMethod::gradient);
    t = markOutlier(t, 10);
    t = dilate(t, 1, 3);
    t = erode(t, 1, 3);
    t = dilate(t, 3, 1);
    t = erode(t, 3, 1);
    std::cout << "get contour.." << std::endl;
    auto ttmp = getContours(t, eg::contourMethod::suzuki);
    Mat2d out(png.info.height, png.info.width);
    out.setConstant(0);
    for(int i = 0; i < ttmp.first.size(); i++) {
        Path tmp = ttmp.first[i];
        for(int j = 1; j < tmp.size(); j++) {
            Segment tttmp = std::make_pair(tmp[j - 1], tmp[j]);
            drawSegmentDirect(out, tttmp, 255);
        }
    }
    i = mat2dToImage(out);
    png.setImage(i);
    png.saveImage("approx-" + inputPath);
}
