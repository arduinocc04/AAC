/**
 * @file main.cpp
 * @author Daniel cho
 * @date 2023.11.18
 * @version 0.0.1
 */
#include <iostream>

#include "egLoader.hpp"
#include "egProcessing.hpp"

using namespace eg::imgproc;

#define ASCII_WIDTH 16
#define ASCII_HEIGHT 22

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use Program proeprly." << std::endl;
    }

    std::string inputPath = argv[1];
    eg::PNG png;
    png.openImage(inputPath);

    png.divideImageByLength(ASCII_HEIGHT, ASCII_WIDTH);
    int gcCnt = png.getGridColCnt();
    int grCnt = png.getGridRowCnt();

    std::cout << grCnt << "x" << gcCnt << std::endl;

    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            Image raw = png.getImageAtGrid(i, j);
            Mat2d sample = cvtGray(raw, eg::grayCvtMethod::mean);
            Eigen::Tensor<double, 0> t = sample.sum();
            double ratio = t(0)/(sample.dimensions()[0]*sample.dimensions()[1]);
            if(ratio < 30) std::cout << " ";
            else if(ratio < 70) std::cout << ".";
            else if(ratio < 110) std::cout << "^";
            else if(ratio < 150) std::cout << "o";
            else if(ratio < 190) std::cout << "*";
            else std::cout << "@";
        }
        std::cout << std::endl;
    }
}
