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

std::vector<Path> getContours(Mat2d & bin, int method);

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
}
