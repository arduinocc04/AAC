/**
 * @file PNGMath.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "PNGMath.hpp"

using namespace Eigen;
using Mat2d = Tensor<double, 2>;

Mat2d eg::math::conv(Mat2d & input, Mat2d & kernel) {
    int kh = kernel.dimensions()[0];
    int kw = kernel.dimensions()[1];
    int h = input.dimensions()[0];
    int w = input.dimensions()[1];
    if((kh % 2 == 0) || (kw % 2 == 0)) // shape of kernel must odd * odd.
        throw eg::exceptions::InvalidParameter();

    Mat2d out(h, w);
    out.setConstant(0);

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            for(int dx = -kh/2; dx < kh/2 + 1; dx++) {
                if(i + dx < 0 || i + dx >= h)
                    continue;
                for(int dy = -kw/2; dy < kw/2 + 1; dy++) {
                    if(j + dy < 0 || j + dy >= w)
                        continue;
                    out(i, j) += kernel(dx + kh/2, dy + kw/2)*input(i + dx, j + dy);
                }
            }
        }
    }
    return out;
}

double eg::math::rmse(Mat2d & a, Mat2d & b) {
    if(a.dimensions()[0] != b.dimensions()[0] ||
        a.dimensions()[1] != b.dimensions()[1])
            throw exceptions::InvalidParameter();
    Tensor<double, 0> tmp = (a-b).square().sum();
    return std::sqrt(tmp(0));
}

Mat2d eg::math::inflate(Mat2d & a, int h, int w) {
    int ah = a.dimensions()[0], aw = a.dimensions()[1];
    if(ah > h || aw > w)
        throw exceptions::InvalidParameter();
    if(ah == h && aw == w)
        return a;

    Mat2d res(h, w);
    res.setConstant(0);
    for(int i = 0; i < ah; i++) {
        for(int j = 0; j < aw; j++) {
                res(i, j) = a(i, j);
        }
    }

    return res;
}
