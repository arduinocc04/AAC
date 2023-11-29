/**
 * @file get-contours.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include <iostream>
#include "egProcessing.hpp"
#include "egLoader.hpp"

using namespace eg::imgproc;

void print(eg::Mat2d & a) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            std::cout << a(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

void print(std::pair<Paths, std::vector<int>> & a) {
    for(int i = 0; i < a.first.size(); i++) {
        for(int j = 0; j < a.first[i].size(); j++) {
            std::cout << a.first[i][j].first << " " << a.first[i][j].second << std::endl;
        }
    }
}

int main(int argc, char * argv[]) {
    if(argc != 2) {
        eg::Mat2d t(10, 10);
        t.setConstant(0);
        for(int i = 1; i < 8; i++) t(1, i) = 1;
        for(int i = 1; i < 6; i++) t(i, 1) = 1;
        for(int i = 1; i < 6; i++) t(i, 3) = 1;
        for(int i = 1; i < 6; i++) t(i, 7) = 1;
        for(int i = 1; i < 8; i++) t(6, i) = 1;
        t(3, 9) = 1;
        print(t);

        std::pair<Paths, std::vector<int>> ans = getContours(t, eg::contourMethod::suzuki);
        std::cout << std::endl << "======ans=====\n";
        print(ans);
        return 0;
    }

    PNG png;
    std::string inputPath = argv[1];
    png.openImage(inputPath);
    Image i = png.copy();
    std::cout << "Cvt gray Image..." << std::endl;
    Mat2d t = cvtGray(i, grayCvtMethod::mean);
    std::cout << "Get edge Image..." << std::endl;
    t = getEdge(t, edgeDetectMethod::gradient);
    t = markOutlier(t, 10);
    t = dilate(t, 1, 3);
    t = erode(t, 1, 3);
    t = dilate(t, 3, 1);
    t = erode(t, 3, 1);
    /*
    for(int i = 0; i < 10; i++) {
        std::cout << i + 1 << "/10 blurring Image..." << std::endl;
        t = blur(t, blurMethod::gaussian);
    }
    */
    std::cout << "get contour.." << std::endl;
    auto ttmp = getContours(t, eg::contourMethod::suzuki);
    Mat2d ans(png.info.height, png.info.width);
    ans.setConstant(0);
    for(int i = 0; i < ttmp.first.size(); i++) {
        for(int j = 0; j < ttmp.first[i].size(); j++) {
            ans(ttmp.first[i][j].first, ttmp.first[i][j].second) = 255;
        }
    }
    i = mat2dToImage(ans);
    png.setImage(i);
    png.saveImage("contour-" + inputPath);
}
