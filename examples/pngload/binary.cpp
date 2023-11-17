#include <iostream>
#include "PNGLoader.hpp"

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }
    eg::PNG png;
    png.openImage(argv[1]);
    png.cvtGray(eg::grayCvtMethod::mean);
    png.getEdge(eg::edgeDetectMethod::gradient);
    for(int i = 0; i < 10; i++) {
        std::cout << "\r" << i + 1 << "/10 blurring input image";
        png.blur(eg::blurMethod::gaussian);
    }
    png.binary(10);
    png.saveImage("bin-" + std::string(argv[1]));
}
