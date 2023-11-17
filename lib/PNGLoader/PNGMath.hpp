#ifndef __PNGEXCEPTIONS_H
#include "PNGExceptions.hpp"
#endif
#define __PNGEXCEPTIONS_H

#ifndef __Tensor
#include <unsupported/Eigen/CXX11/Tensor>
#endif
#define __Tensor

namespace eg::math {
/**
 * @todo Study other's implementation. I don't know how to behave at corners.
 */
Eigen::Tensor<double, 2> conv(
                    Eigen::Tensor<double, 2> & input,
                    Eigen::Tensor<double, 2> & kernel
                    );

/**
 * @brief calculate rmse of same size 2d tensor
 * @attention two tensor must have same size
 * @param a 2d double tensor
 * @param b 2d double tensor
 * @return rmse double type
 */
double rmse(Eigen::Tensor<double, 2> & a, Eigen::Tensor<double, 2> & b);

/**
 * @brief size up tensor by fill zeros
 * @param a 2d double tensor
 * @param h target height
 * @param w target width
 * @return scaled tensor
 */
Eigen::Tensor<double, 2> inflate(Eigen::Tensor<double, 2> & a, int h, int w);

}
