#include "PNGMath.hpp"

Eigen::Tensor<double, 2> eg::math::conv(Eigen::Tensor<double, 2> & input, Eigen::Tensor<double, 2> & kernel) {
    int kx = kernel.dimensions()[0];
    int ky = kernel.dimensions()[1];
    int h = input.dimensions()[0];
    int w = input.dimensions()[1];
    if((kx%2)*(ky%2) == 0) // shape of kernel must odd * odd.
        throw eg::exceptions::InvalidParameter();

    Eigen::Tensor<double, 2> out(h, w);
    out.setConstant(0);

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            for(int dx = -kx/2; dx < kx/2 + 1; dx++) {
                if(i + dx < 0 || i + dx >= h)
                    continue;
                for(int dy = -ky/2; dy < ky/2 + 1; dy++) {
                    if(j + dy < 0 || j + dy >= w)
                        continue;
                    out(i, j) += kernel(dx + kx/2, dy + ky/2)*input(i + dx, j + dy);
                }
            }
        }
    }
    return out;
}
