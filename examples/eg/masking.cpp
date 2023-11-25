#include <iostream>

#include "egLoader.hpp"
#include "egProcessing.hpp"

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
    if(argc != 2)
        std::cout << "Use Program properly!" << std::endl;

    eg::PNG png;
    std::string inputPath = argv[1];
    png.openImage(inputPath);
    Image i = png.copy();
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    t = getEdge(t, eg::edgeDetectMethod::gradient);
    t = 255*getMask(t);
    i = mat2dToImage(t);
    png.setImage(i);
    png.saveImage("mask-"+inputPath);
}
