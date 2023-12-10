/**
 * @file egMath.hpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#ifndef __EGEXCEPTIONS_H
#include "egExceptions.hpp"
#endif
#define __EGEXCEPTIONS_H

#ifndef __Tensor
#include <unsupported/Eigen/CXX11/Tensor>
#endif
#define __Tensor

#ifndef __EGTYPES_H
#include "egTypes.hpp"
#endif
#define __EGTYPES_H

#include "egProcessing.hpp"

using namespace eg;

#include <vector>

namespace eg::math {
/**
 * @todo Study other's implementation. I don't know how to behave at corners.
 */
Eigen::Tensor<double, 2> conv(const Mat2d & input, const Mat2d & kernel);

double compareMat2d(const Mat2d & a, const Mat2d & b, int method);

double compareHistogram(const Mat2d & ha, const Mat2d & hb, int method);

double calcDeformLocal(const Segment & before, const Segment & after);

double calcAccess(const Segment & before, const Segment & after, const Segments & ss, const Segments & original, const std::vector<int> & candidates);

double calcDeform(const Segment & before, const Segment & after, const Segments & ss, const Segments & original, const std::vector<int> & candidates);
}

