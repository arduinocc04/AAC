/**
 * @file egProcessing.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "egProcessing.hpp"
#include "egTypes.hpp"
#include "egGeometry.hpp"
#include "egMath.hpp"
#include <cmath>

using namespace eg;
using namespace eg::imgproc;
using namespace eg::geo;

Mat2d eg::imgproc::binary(Mat2d & gray, double threshold) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Mat2d bin(h, w);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            bin(i, j) = (double)(gray(i, j) > threshold);

    return bin;
}

Mat2d eg::imgproc::getMask(Mat2d & gray) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Dots candidates;
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            if(gray(i, j))
                candidates.push_back(std::make_pair(i, j));
    Dots hull = getConvexHull(candidates);

    Mat2d ans = Mat2d(h, w);
    ans.setConstant(0);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(isPointInsideHull(std::make_pair(i, j), hull, true))
                ans(i, j) = 1;
        }
    }
    return ans;
}

Mat2d blurGaussian(Mat2d & gray) {
    Mat2d kernel(3, 3);
    kernel(0, 0) = (double)1/16, kernel(0, 1) = (double)2/16, kernel(0, 2) = (double)1/16;
    kernel(1, 0) = (double)2/16, kernel(1, 1) = (double)4/16, kernel(1, 2) = (double)2/16;
    kernel(2, 0) = (double)1/16, kernel(2, 1) = (double)2/16, kernel(2, 2) = (double)1/16;
    return eg::math::conv(gray, kernel);
}

Mat2d eg::imgproc::blur(Mat2d & gray, int method) {
    switch(method) {
        case eg::blurMethod::gaussian:
            return blurGaussian(gray);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d cvtGrayMean(Image & imageRGBA) {
    if(!imageRGBA.data())
        throw eg::exceptions::ImageNotOpened();

    int h = imageRGBA.dimensions()[0];
    int w = imageRGBA.dimensions()[1];
    Mat2d gray(h, w);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            double mean = 0;
            for(int k = 0; k < 3; k++)
                mean += imageRGBA(i, j, k);
            mean /= 3;
            gray(i, j) = mean;
        }
    }
    return gray;
}

Mat2d eg::imgproc::cvtGray(Image & imageRGBA, int method) {
    switch(method) {
        case eg::grayCvtMethod::mean:
            return cvtGrayMean(imageRGBA);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d getEdgeGrad(Mat2d & gray) {
    Eigen::Tensor<double, 2> kernel(3, 3);
    kernel(0, 0) = -1, kernel(0, 1) = -1, kernel(0, 2) = -1;
    kernel(1, 0) = -1, kernel(1, 1) =  8, kernel(1, 2) = -1;
    kernel(2, 0) = -1, kernel(2, 1) = -1, kernel(2, 2) = -1;

    return eg::math::conv(gray, kernel);
}

Mat2d eg::imgproc::getEdge(Mat2d & gray, int method) {
    switch(method) {
        case eg::edgeDetectMethod::gradient:
            return getEdgeGrad(gray);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d extCenterlineGrassfire(Mat2d & bin) {
    Mat2d mask = imgproc::getMask(bin);
    return eg::math::grassfire(bin, mask);
}

Mat2d eg::imgproc::extractCenterline(Mat2d & bin, int method) {
    switch(method) {
        case eg::centerlineMethod::grassfire:
            return extCenterlineGrassfire(bin);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d addZeroPadding(Mat2d & a, int h, int j, int k, int l) {
    int ah = a.dimensions()[0];
    int aw = a.dimensions()[1];
    Mat2d res(ah + j + k, aw + h + l);
    res.setConstant(0);
    for(int r = 0; r < ah; r++)
        for(int c = h; c < aw; c++)
            res(r + k, c + h) = a(r, c);
    return res;
}

Mat2d eg::imgproc::addPadding(Mat2d & a, int h, int j, int k, int l, int method) {
    if(h < 0 || j < 0 || k < 0 || l < 0)
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::paddingMethod::zero:
            return addZeroPadding(a, h, j, k, l);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d eg::imgproc::inflate(Mat2d & a, int h, int w) {
    int ah = a.dimensions()[0], aw = a.dimensions()[1];
    return imgproc::addPadding(a, 0, h - ah, 0, w - aw,
                     eg::paddingMethod::zero);
}

unsigned char saturate(double x) {
    if(x < 0) return 0;
    if(x > 255) return 255;
    return (unsigned char)x;
}

Image eg::imgproc::mat2dToImage(Mat2d & a) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    Image res(h, w, 4);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++) {
            for(int k = 0; k < 3; k++)
                res(i, j, k) = std::round(a(i, j));
            res(i, j, 3) = 255;
        }
    return res;
}
