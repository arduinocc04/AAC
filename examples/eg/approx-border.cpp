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
    /*for(int i = 0; i < 10; i++) {
        std::cout << i + 1 << "/10 blurring Image..." << std::endl;
        t = blur(t, blurMethod::gaussian);
    }*/
    t = markOutlier(t, 10);
    /*
    for(int i = 0; i < 10; i++) {
        std::cout << i + 1 << "/10 image opening" << std::endl;
        t = erode(t, 1, 3);
        t = dilate(t, 1, 3);
    }
    for(int i = 0; i < 10; i++) {
        std::cout << i + 1 << "/10 image opening" << std::endl;
        t = erode(t, 3, 1);
        t = dilate(t, 3, 1);
    }
    t = dilate(t, 3, 1);
    t = dilate(t, 1, 3);
    */
    t = dilate(t, 1, 3);
    t = erode(t, 1, 3);
    t = dilate(t, 3, 1);
    t = erode(t, 3, 1);
    std::cout << "get contour.." << std::endl;
    Mat2d out(png.info.height, png.info.width);
    out.setConstant(0);
    auto ttmp = getContours(t, eg::contourMethod::suzuki);
    for(int i = 0; i < ttmp.first.size(); i++) {
        if(!ttmp.first[i].size()) continue;
        std::cout << "started decomposing" << std::endl;
        Segments tmp = eg::trace::decomposePathToSegments(ttmp.first[i], eg::pathDecomMethod::greedy);
        std::cout << "decomposing finished." << std::endl;
        std::cout << "start drawing " << tmp.size() << " segments" << std::endl;
        for(int j = 0; j < tmp.size(); j++) {
            drawSegmentDirect(out, tmp[j], 255);
        }
        std::cout << "end drawing" << std::endl;
    }
    i = mat2dToImage(out);
    png.setImage(i);
    png.saveImage("approx-" + inputPath);
}
