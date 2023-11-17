#include "PNGLoader.hpp"

void eg::PNG::binary(double threshold) {
    if(!info.initialized || !image.data())
        throw exceptions::ImageNotOpened();

    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            playground(i, j) = (playground(i, j) > threshold)*255;
        }
    }
    copyPlaygroundToImage();
}
