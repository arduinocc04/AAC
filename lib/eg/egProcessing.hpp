/**
 * @file egProcessing.hpp
 * @autor Daniel Cho
 * @date 2023.11.24
 * @version 0.0.1
 */
#include <vector>

#include "unsupported/Eigen/CXX11/Tensor"

#include "egMethods.hpp"
#include "egTypes.hpp"

#ifndef __EGEXCEPTIONS_H
#include "egExceptions.hpp"
#endif
#define __EGEXCEPTIONS_H

using namespace eg;

namespace eg::imgproc {
 /**
 * @brief convert opened image to grayscale.
 * @attention you need open image before call this function. Call this method will change playground and opened image.
 * @param method an integer
 * @see cvtGray.cpp
 * @see eg::grayCvtMethod
 * @throws GetMetadataFailed
 * @throws InvalidFormat
 * @throws InvalidParameter
 */
Mat2d cvtGray(Image & imageRGBA, int method);

 /**
  * @brief get Edge of playground.
  * @attention This will change playground and opened image.
  */
Mat2d getEdge(Mat2d & gray, int method);

/**
 * @brief blur playground.
 * @attention This will change playground and opened image
 */
Mat2d blur(Mat2d & gray, int method);

Mat2d binary(Mat2d & gray, double threshold);

Mat2d extractCenterline(Mat2d & bin, int method);

Mat2d saturate(Mat2d & gray);

Mat2d markOutlier(Mat2d & gray, double threshold);

/**
 * @attention param must be binary image.
 */
Mat2d reverse(Mat2d & bin);

Mat2d erode(Mat2d & bin, int kh, int kw);

Mat2d dilate(Mat2d & bin, int kh, int kw);

/**
 * @breif getcontours using algorithm made by Suzuki and Abe
 * @see Topological Structural Analysis of Digitized Binary Images by Border Following(1985) Appendix 1.
 * @param bin binary image.
 * @todo refactor
 */
std::pair<Paths, std::vector<int>> getContours(Mat2d & bin, int method);

Mat2d getMask(Mat2d & gray);

Mat2d addPadding(Mat2d & a, int h, int j, int k, int l, int method);

/**
 * @brief size up tensor by fill zeros
 * @param a 2d double tensor
 * @param h target height
 * @param w target width
 * @return scaled tensor
 */
Mat2d inflate(Mat2d & a, int h, int w);

Image mat2dToImage(Mat2d & a);

Mat2d logpolar(Dots & a, Dot & p);

Mat2d logpolarAll(Dots & a);
/**
 * @attention the border of mask must zero. If not, if will raise segfault.
 */
Mat2d grassfire(Mat2d & a, Mat2d & mask);

Mat2d cycle(const Mat2d & a, int stride);

}
