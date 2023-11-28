/**
 * @file egMath.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include <cmath>

#include "egMath.hpp"
#include "egTypes.hpp"
#include "egProcessing.hpp"

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

// I wanted to name this function as rmse, but stupid compiler thought it's ambiguous.
double calcrmse(Mat2d & a, Mat2d & b) {
    Tensor<double, 0> tmp = (a-b).square().sum();
    return std::sqrt(tmp(0));
}

double bhattacharyyaDist(Mat2d & ha, Mat2d & hb) {
    int h = ha.dimensions()[0];
    int w = ha.dimensions()[1];
    Eigen::Tensor<double, 0> sa = ha.sum();
    Eigen::Tensor<double, 0> sb = hb.sum();
    double n = std::sqrt(sa(0))*std::sqrt(sb(0));
    double s = 0;
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            s += std::sqrt(ha(i, j))*std::sqrt(hb(i, j));
        }
    }
    return std::sqrt(std::abs(1 - 1/n*s));
}

double eg::math::compareHistogram(Mat2d & ha, Mat2d & hb, int method) {
    if(ha.dimensions() != hb.dimensions())
        throw eg::exceptions::InvalidParameter();

    switch(method) {
        case eg::histCmpMethod::bhattacharyya:
            return bhattacharyyaDist(ha, hb);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Dots merge(Paths & a) {
    Dots p;
    for(int i = 0; i < a.size(); i++)
        for(int j = 0; j < a[i].size(); j++)
            p.push_back(a[i][j]);
    return p;
}

double eg::math::compareMat2d(Mat2d & a, Mat2d & b, int method) {
    if(a.dimensions() != b.dimensions())
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::matCmpMethod::rmse:
            return calcrmse(a, b);
        case eg::matCmpMethod::logpolar: {
            std::pair<Paths, std::vector<int>> tmpA = eg::imgproc::getContours(a, eg::contourMethod::suzuki);
            std::pair<Paths, std::vector<int>> tmpB = eg::imgproc::getContours(b, eg::contourMethod::suzuki);
            Dots dotsConsistsA = merge(tmpA.first);
            Dots dotsConsistsB = merge(tmpB.first);
            Mat2d histogramA = eg::imgproc::logpolar(dotsConsistA);
            Mat2d histogramB = eg::imgproc::logpolar(dotsConsistB);
            return compareHistogram(histogramA, histogramB, eg::histCmpMethod::bhattacharyya);
        }
        default:
            throw eg::exceptions::InvalidParameter();
    }
}
