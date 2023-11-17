#include "PNGLoader.hpp"

void eg::PNG::getEdgeGrad() {
    Eigen::Tensor<double, 2> kernel(3, 3);
    kernel(0, 0) = -1, kernel(0, 1) = -1, kernel(0, 2) = -1;
    kernel(1, 0) = -1, kernel(1, 1) =  8, kernel(1, 2) = -1;
    kernel(2, 0) = -1, kernel(2, 1) = -1, kernel(2, 2) = -1;

    playground = math::conv(playground, kernel);
}

void eg::PNG::getEdge(int method) {
    if(!info.initialized)
        throw exceptions::GetMetadataFailed();
    if(!playground.data())
        throw exceptions::InvalidFormat();

    switch(method) {
        case edgeDetectMethod::gradient:
            getEdgeGrad();
            break;
        default:
            throw exceptions::InvalidParameter();
            break;
    }

    copyPlaygroundToImage();
}
