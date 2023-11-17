#include <iostream>
#include "PNGLoader.hpp"

int main(int argc, char * argv[]) {
    if(argc != 3) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }
    eg::PNG png;
    png.openImage(argv[1]);
    png.cvtGray(eg::grayCvtMethod::mean);
    png.getEdge(eg::edgeDetectMethod::gradient);
    int n = std::stoi(argv[2]);
    for(int i = 0; i < n; i++) png.blur(eg::blurMethod::gaussian);
    png.saveImage("ed+blur" + std::to_string(n) + "-" + std::string(argv[1]));
}
