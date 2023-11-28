/**
 * @file histogram.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include <iostream>
#include "egLoader.hpp"
#include "egProcessing.hpp"
#include "egMath.hpp"

using namespace eg;
using namespace eg::imgproc;

void print(Mat2d & a) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            std::cout << a(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char * argv[]) {
    PNG png;
    png.openImage(argv[1]);
    Image i = png.copy();
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    // t = getEdge(t, eg::edgeDetectMethod::gradient);
    t = markOutlier(t, 70);
    std::cout << "=========Image 1" << std::endl;
    print(t);
    std::cout << "==========" << std::endl;
    std::pair<Paths, std::vector<int>> tmp = getContours(t, eg::contourMethod::suzuki);
    Path dots1;
    for(int i = 0; i < tmp.first.size(); i++) {
        for(int j = 0; j < tmp.first[i].size(); j++)
            dots1.push_back(tmp.first[i][j]);
    }

    std::cout << "Generaing logpolar of dots1: " << dots1.size() << std::endl;
    Mat2d h1 = eg::imgproc::logpolar(dots1);
    std::cout << "generated\n";

    png.openImage(argv[2]);
    i = png.copy();
    t = cvtGray(i, eg::grayCvtMethod::mean);
    //t = getEdge(t, eg::edgeDetectMethod::gradient);
    t = markOutlier(t, 70);
    std::cout << "============Image 2" << std::endl;
    print(t);
    std::cout << "============" << std::endl;
    tmp = getContours(t, eg::contourMethod::suzuki);
    Path dots2;
    for(int i = 0; i < tmp.first.size(); i++)
        for(int j = 0; j < tmp.first[i].size(); j++)
            dots2.push_back(tmp.first[i][j]);

    std::cout << "Generaing logpolar of dots2: " << dots2.size() << std::endl;
    Mat2d h2 = eg::imgproc::logpolar(dots2);
    std::cout << "generated\n";

    std::cout << "===========Histo 1" << std::endl;
    print(h1);
    std::cout << "===========Histo 2" << std::endl;
    print(h2);
    std::cout << "==========" << std::endl;
    double err = eg::math::compareHistogram(h1, h2, eg::histCmpMethod::bhattacharyya);
    std::cout << err << std::endl;
}
