/**
 * @file egProcessing.hpp
 * @autor Daniel Cho
 * @date 2023.11.24
 * @version 0.0.1
 */
#include <vector>

#include "unsupported/Eigen/CXX11/Tensor"

#ifndef __EGMETHODS_H
#include "egMethods.hpp"
#endif
#define __EGMETHODS_H

#ifndef __EGTYPES_H
#include "egTypes.hpp"
#endif
#define __EGTYPES_H

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
Mat2d cvtGray(const Image & imageRGBA, int method);

 /**
  * @brief get Edge of playground.
  * @attention This will change playground and opened image.
  */
Mat2d getEdge(const Mat2d & gray, int method);

/**
 * @brief blur playground.
 * @attention This will change playground and opened image
 */
Mat2d blur(const Mat2d & gray, int method);

Mat2d binary(const Mat2d & gray, double threshold);

Mat2d extractCenterline(const Mat2d & bin, int method);

Mat2d saturate(const Mat2d & gray);

/**
 * @attention This method also handles negative-valued pixels.
 */
Mat2d markOutlier(const Mat2d & gray, double threshold);

/**
 * @attention param must be binary image.
 */
Mat2d reverse(const Mat2d & bin);

Mat2d erode(const Mat2d & bin, int kh, int kw);

Mat2d dilate(const Mat2d & bin, int kh, int kw);

/**
 * @breif getcontours using algorithm made by Suzuki and Abe
 * @see Topological Structural Analysis of Digitized Binary Images by Border Following(1985) Appendix 1.
 * @param bin binary image.
 * @todo refactor
 */
std::pair<Paths, std::vector<int>> getContours(const Mat2d & bin, int method);

Mat2d getMask(const Mat2d & gray);

Mat2d addPadding(const Mat2d & a, int h, int j, int k, int l, int method);

/**
 * @brief size up tensor by fill zeros
 * @param a 2d double tensor
 * @param h target height
 * @param w target width
 * @return scaled tensor
 */
Mat2d inflate(const Mat2d & a, int h, int w);

Image mat2dToImage(const Mat2d & a);

Mat2d logpolar(const Dots & a, const Dot & p);

Mat2d logpolarAll(const Dots & a);

/**
 * @attention the border of mask must zero. If not, if will raise segfault.
 */
Mat2d grassfire(const Mat2d & a, const Mat2d & mask);

Mat2d cycle(const Mat2d & a, int stride);

/**
 * @breif draw segment using Bresenham's line algorithm
 * @see https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
 */
Mat2d drawSegment(const Mat2d & a, const Segment & s, int val);

void drawSegmentDirect(Mat2d & a, const Segment & s, int val);

Mat2d drawSegments(const Mat2d & a, const Segments & ss, int val);

/**
 * @return return binary image. (0 or 1)
 */
Mat2d approxUsingSegments(const Mat2d & a);

Mat2d resizeImage(const Mat2d & bin, int method, int targetH, int targetW);
}
