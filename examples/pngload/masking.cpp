#include <iostream>

#include "PNGLoader.hpp"

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
    png.cvtGray(eg::grayCvtMethod::mean);
    png.getEdge(eg::edgeDetectMethod::gradient);
    Eigen::Tensor<double, 2> mask = eg::math::getMask(*png.getPlayground());
    for(int i = 0; i < mask.dimensions()[0]; i++) {
        for(int j = 0; j < mask.dimensions()[1]; j++) {
            (*png.getPlayground())(i, j) = mask(i, j)*255;
        }
    }
    png.binary(150); // called this to copy playground to image.
    png.saveImage("mask-"+inputPath);
}
