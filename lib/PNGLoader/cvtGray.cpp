/**
 * @file cvtGray.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "PNGLoader.hpp"

void eg::PNG::cvtGrayMean() {
    if(!info.initialized || !image.data())
        throw exceptions::ImageNotOpened();

    for(int i = 0; i < info.height; i++) {
        for(int j = 0; j < info.width; j++) {
            int mean = 0;
            for(int k = 0; k < 3; k++)
                mean += image(i, j, k);
            mean /= 3;
            image(i, j, 0) = image(i, j, 1) = image(i, j, 2) = mean;
        }
    }
}

void eg::PNG::cvtGray(int method) {
    if(!info.initialized)
        throw exceptions::GetMetadataFailed();
    if(info.colorType != PNG_COLOR_TYPE_RGB_ALPHA)
        throw exceptions::InvalidFormat();
    switch(method) {
        case grayCvtMethod::mean:
            cvtGrayMean();
            break;
        default:
            throw exceptions::InvalidParameter();
            break;
    }
    allocPlayground();
    copyImageToPlayground();
}
