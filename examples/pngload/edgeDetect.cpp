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
    png.saveImage("edge-" + std::string(argv[1]));
}