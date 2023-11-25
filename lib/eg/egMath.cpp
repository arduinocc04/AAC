/**
 * @file PNGMath.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "egMath.hpp"
#include "egTypes.hpp"

using namespace eg;
using namespace Eigen;

/**
 * @todo padding 붙여서 구현하는게 더 깔끔할듯
 */
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

/**
 * @attention the border of mask must zero. If not, if will raise segfault.
 */
Mat2d eg::math::grassfire(Mat2d & a, Mat2d & mask) {
    if(a.dimensions() != mask.dimensions())
        throw exceptions::InvalidParameter();
    Mat2d ans(a.dimensions());
    ans.setConstant(0);
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(mask(i, j) && a(i, j) > 0)
                ans(i, j) = 1 + std::min(ans(i - 1, j), ans(i, j - 1)); // use taxi distance
        }
    }
    for(int i = h - 1; i >= 0; i--) {
        for(int j = w - 1; j >= 0; j--) {
            if(mask(i, j) && a(i, j) > 0)
                ans(i, j) = std::min(ans(i, j),
                                     1 + std::min(ans(i + 1, j), ans(i, j + 1)));
        }
    }

    return ans;
}

