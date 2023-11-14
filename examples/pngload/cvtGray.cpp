#include <iostream>
#include "PNGLoader.hpp"

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }
    eg::PNG * png = new eg::PNG;
    png->openImage(argv[1]);
    png->cvtGray(eg::grayCvtMethod::mean);
    png->saveImage("gray-" + std::string(argv[1]));
    delete png;
}
