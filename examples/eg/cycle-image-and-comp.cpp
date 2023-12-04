/**
 * @file cycle-image-and-comp.cpp
 * @author Daniel Cho
 * @date 2023.11.29
 * @version 0,0,1
 */
#include <iostream>

#include "egProcessing.hpp"
#include "egLoader.hpp"
#include "egMath.hpp"

using namespace eg;
using namespace eg::imgproc;

void print(Eigen::Tensor<double, 2> & a) {
    for(int i = 0; i < a.dimensions()[0]; i++) {
        for(int j = 0; j < a.dimensions()[1]; j++) {
            std::cout << (a(i, j) > 0);
        }
        std::cout << std::endl;
    }
}

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Run program properly!" << std::endl;
        return -1;
    }
    PNG png;
    png.openImage(argv[1]);
    Image i = png.copy();
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    t = binary(t, 70);
    int w = t.dimensions()[1];
    std::cout << "Original======\n";
    print(t);
    for(int j = 1; j < w; j++) {
        std::cout << "Cycled=====" << j << "\n";
        Mat2d tmp = cycle(t, j);
        print(tmp);
        double err = eg::math::compareMat2d(t, tmp, eg::matCmpMethod::logpolar);
        std::cout << err << std::endl;
    }
}
