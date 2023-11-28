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

#include "egTypes.hpp"
using namespace eg;

#include <vector>

namespace eg::math {
/**
 * @todo Study other's implementation. I don't know how to behave at corners.
 */
Eigen::Tensor<double, 2> conv(
                    Eigen::Tensor<double, 2> & input,
                    Eigen::Tensor<double, 2> & kernel
                    );

double compareMat2d(Mat2d & a, Mat2d & b, int method);

double compareHistogram(Mat2d & ha, Mat2d & hb, int method);
}

