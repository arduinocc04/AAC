#include "PNGLoader.hpp"

void eg::PNG::blur(int method) {
    if(!info.initialized || !image.data())
        throw exceptions::ImageNotOpened();
    if(!playground.data())
        throw exceptions::InvalidFormat();

    switch(method) {
        case blurMethod::gaussian:
            blurGaussian();
            break;
        default:
            throw exceptions::InvalidParameter();
    }
    copyPlaygroundToImage();
}

void eg::PNG::blurGaussian() {
    Eigen::Tensor<double, 2> kernel(3, 3);
    kernel(0, 0) = (double)1/16, kernel(0, 1) = (double)2/16, kernel(0, 2) = (double)1/16;
    kernel(1, 0) = (double)2/16, kernel(1, 1) = (double)4/16, kernel(1, 2) = (double)2/16;
    kernel(2, 0) = (double)1/16, kernel(2, 1) = (double)2/16, kernel(2, 2) = (double)1/16;
    playground = math::conv(playground, kernel);
}
